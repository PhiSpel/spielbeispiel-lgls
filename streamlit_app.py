# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:38:24 2022

@author: Philipp Spelten
"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcol
import matplotlib.cm as cm

import plotly.express as px
#import chart_studio.plotly as py
import plotly.graph_objects as go
import plotly
#import plotly.colors.sequential as seq_color
from plotly.figure_factory import create_quiver
# Identical to Adam's answer
import plotly.colors
from PIL import ImageColor

import streamlit as st
from streamlit_plotly_events import plotly_events

import math

###############################################################################
# CALCULATIONS
###############################################################################

# maybe use @st.cache()?
def update_data(nodes,members,support,f_ext):

    n_nodes = len(nodes)
    nodes_range = np.arange(0,n_nodes)
    n_members = len(members)
    members_range = np.arange(0,n_members)
    nsupport = len(support)
    support_range = np.arange(0,nsupport)
    nf_ext = len(f_ext)
    force_range = np.arange(0,nf_ext)
    
    members = np.append(members, np.zeros([n_members,3]),axis = 1)
    #members[:,0:2] = members[:,0:2].astype(int)
    # compute angles
    for m in members_range:
        node1 = int(members[m,0]-1)
        node2 = int(members[m,1]-1)
        stab = nodes[node2,:]-nodes[node1,:]
        direction = np.sign(stab[1])
        if direction == 0:
            direction = 1
        alpha = np.arccos(stab[0]/np.linalg.norm(stab))*direction
        members[m,2] = np.cos(alpha)
        members[m,3] = np.sin(alpha)
        members[m,4] = alpha
    
    # compute matrix coefficients
    ks1 = np.zeros([n_nodes,n_members])
    ks2 = np.zeros([n_nodes,n_members])
    matrix = np.zeros([2*n_nodes,n_members])
    for k in nodes_range:
        # compute x and y forces per member
        for m in members_range:
            node1 = int(members[m,0]-1)
            node2 = int(members[m,1]-1)
            ks1[k,m] = int(k==node1)
            ks2[k,m] = -int(k==node2)
        ks = ks1+ks2
        # forces along x
        matrix[k,:] = np.multiply(ks[k,:],members[:,2])
        # forces along y
        matrix[k+n_nodes,:] = np.multiply(ks[k,:],members[:,3])
        
    # add support forces to matrix
    matrix = np.append(matrix, np.zeros([len(matrix),3]),axis = 1)
    for s in support_range:
        node1 = int(support[s,0]-1)
        angle = support[s,1]
        matrix[node1,n_members+s] = np.cos(angle)
        matrix[node1+n_nodes,n_members+s] = np.sin(angle)
        
    # compute right hand side
    rhs = np.zeros([2*n_nodes,1])
    for force in force_range:
        node_index = int(f_ext[force,0]-1)
        angle = f_ext[force,1]
        newtons = f_ext[force,2]
        # external force along x
        rhs[node_index,0] = -newtons*np.cos(angle)
        # external force along y
        rhs[node_index+n_nodes,0] = -newtons*np.sin(angle)
        
    # solve for unknowns
    internal_forces = np.linalg.solve(matrix, rhs)
    
    return matrix, rhs, internal_forces

def new_member(new_nodes):
    node1 = new_nodes[0]
    node2 = new_nodes[1]
    st.session_state.members = np.append(st.session_state.members,[[node1,node2]],axis=0)
    return

###############################################################################
# PLOTS
###############################################################################

def get_continuous_color(colorscale, intermed):
    """
    Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    color for any value in that range.

    Plotly doesn't make the colorscales directly accessible in a common format.
    Some are ready to use:
    
        colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    Others are just swatches that need to be constructed into a colorscale:

        viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
        colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    :param intermed: value in the range [0, 1]
    :return: color in rgb string format
    :rtype: str
    """
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    hex_to_rgb = lambda c: "rgb" + str(ImageColor.getcolor(c, "RGB"))

    if intermed <= 0 or len(colorscale) == 1:
        c = colorscale[0][1]
        return c if c[0] != "#" else hex_to_rgb(c)
    if intermed >= 1:
        c = colorscale[-1][1]
        return c if c[0] != "#" else hex_to_rgb(c)

    for cutoff, color in colorscale:
        if intermed > cutoff:
            low_cutoff, low_color = cutoff, color
        else:
            high_cutoff, high_color = cutoff, color
            break

    if (low_color[0] == "#") or (high_color[0] == "#"):
        # some color scale names (such as cividis) returns:
        # [[loc1, "hex1"], [loc2, "hex2"], ...]
        low_color = hex_to_rgb(low_color)
        high_color = hex_to_rgb(high_color)

    return plotly.colors.find_intermediate_color(
        lowcolor=low_color,
        highcolor=high_color,
        intermed=((intermed - low_cutoff) / (high_cutoff - low_cutoff)),
        colortype="rgb",
    )

def update_plot(internal_forces,members,nodes,f_ext,support):

    fig = go.Figure()
    
    # draw beams
    def get_color(colorscale_name, loc):
        from _plotly_utils.basevalidators import ColorscaleValidator
        # first parameter: Name of the property being validated
        # second parameter: a string, doesn't really matter in our use case
        cv = ColorscaleValidator("colorscale", "")
        # colorscale will be a list of lists: [[loc1, "rgb1"], [loc2, "rgb2"], ...] 
        colorscale = cv.validate_coerce(colorscale_name)
        
        if hasattr(loc, "__iter__"):
            return [get_continuous_color(colorscale, x) for x in loc]
        return get_continuous_color(colorscale, loc)
        
    max_force=max(abs(internal_forces))[0]
    normed_force = (internal_forces/max_force + 1)/2
    
    for m in np.arange(0,len(members)):
        node1_index = int(members[m,0]-1)
        node2_index = int(members[m,1]-1)
        node1_coord = nodes[node1_index,:]
        node2_coord = nodes[node2_index,:]
        x=[]
        y=[]
        x.append(node1_coord[0])
        x.append(node2_coord[0])
        y.append(node1_coord[1])
        y.append(node2_coord[1])
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            mode="lines + markers",
            name='member #' + str(m+1),
            showlegend=False,
            line=dict(
                color=get_color('rdbu',normed_force[m][0])
                ),
            marker=dict(
                cmax=max_force,
                cmin=-max_force,
                colorbar=dict(
                    title="Force intensity",
                    orientation='h' # this does nothing! :/
                ),
                colorscale="rdbu"
            )
            ))
    
    # create a heatmap to be added below
    forcemap = px.imshow([np.linspace(-max_force,max_force,100)],
                         color_continuous_scale='rdbu')
    forcemap.layout.coloraxis.showscale = False
        
    # draw nodes last, so they can be selected
    fig.add_trace(go.Scatter(x=nodes[:,0],y=nodes[:,1],
                    mode='markers',
                    text=np.arange(1,len(nodes)+1),
                    name='nodes'
                    ))
    
    # calculate coordinates of forces
    x0 = []
    y0 = []
    fx = []
    fy = []
    for f in np.arange(0,len(f_ext)):
        node_index = int(f_ext[f,0]-1)
        x0.append(nodes[node_index,0])
        y0.append(nodes[node_index,1])
        angle = f_ext[f,1]
        newtons = f_ext[f,2]
        # external force along x
        fx.append(newtons*np.cos(angle))
        # external force along y
        fy.append(newtons*np.sin(angle))
    # draw forces  
    quiver = create_quiver(x0,y0,fx,fy,scale=0.05,line=(dict(color='red')),name='Forces')
    fig.add_traces(data=quiver.data)
        
    # draw supports
    support_length = -5
    x0 = []
    y0 = []
    sx = []
    sy = []
    for s in np.arange(0,len(support)):
        node_index = int(support[s,0]-1)
        x0.append(nodes[node_index,0])
        y0.append(nodes[node_index,1])
        angle = support[s,1]
        # external force along x
        sx.append(support_length*np.cos(angle))
        # external force along y
        sy.append(support_length*np.sin(angle))
    # draw forces  
    quiver = create_quiver(x0,y0,sx,sy,
                           scale=0.05,
                           line=(dict(color='green')),
                           name='Supports')
    fig.add_traces(data=quiver.data)
        
    #draw center points to be able to deselect beams
    centresx=[]
    centresy=[]
    member_names=[]
    for m in np.arange(0,len(members)):
        node1_index = int(members[m,0]-1)
        node2_index = int(members[m,1]-1)
        node1_coord = nodes[node1_index,:]
        node2_coord = nodes[node2_index,:]
        centresx.append((node1_coord[0]+node2_coord[0])/2)
        centresy.append((node1_coord[1]+node2_coord[1])/2)
        member_names.append('deselect member #' + str(m+1) + '; current force: ' + str(round(internal_forces[m][0])))
        
    fig.add_trace(go.Scatter(
        x=centresx,
        y=centresy,
        mode='markers',
        marker=dict(symbol='x',color='orange'),
        name='members for deselection',
        text=member_names
        ))
        
    return fig,forcemap

###############################################################################
# PRINTING OUTPUT
###############################################################################

def bmatrix(a,matrixtype=''):
    # if matrixtype == 'h':
    #     text = r' \begin{matrix} '
    #     text += '\n'
    #     for x in range(len(a)):
    #         text+= str(a[x])
    #         text += r' & '
            
    #     text += r' \end{matrix} '
    #     return text
    if matrixtype == 'b':
        text = r' \begin{bmatrix} '
    elif matrixtype == 'v':
        text = r' \begin{matrix} '
    text += '\n'
    for x in range(len(a)):
        for y in range(len(a[x])):
            text += str(a[x][y])
            text += r' & '
        text = text[:-2]
        text += r'\\'
        text += '\n'
    if matrixtype == 'b':
        text += r' \end{bmatrix} '
    elif matrixtype == 'v':
        text += r' \end{matrix} '
    
    return text

def print_equations(matrix, rhs, internal_forces,n_beams,n_bcs,decimals,textsize):
    label_vector = [[r' \text{ ' + str(math.floor(x/2)+1) + '}'] for x in np.arange(0,len(rhs))]
    matrix = np.round(matrix,decimals)
    matrix[matrix==0] = 0
    rhs = np.round(rhs,decimals)
    rhs[rhs==0] = 0
    internal_forces = np.round(internal_forces,decimals)
    internal_forces[internal_forces==0] = 0
    equation_string = r'$'
    equation_string += textsize
    equation_string += r' \begin{matrix} '
    equation_string += r'\text{nodes} & \text{' + str(n_beams) + ' beams and ' + str(n_bcs) + ' boundary conditions}'
    #equation_string += bmatrix(label_vector,'h')
    equation_string += r' & \text{internal forces} & \text{external forces}\\'
    equation_string += bmatrix(label_vector,'v')
    equation_string += ' & '
    equation_string += bmatrix(matrix,'b')
    equation_string += ' & '
    equation_string += '\cdot '
    equation_string += bmatrix(internal_forces,'b')
    equation_string += ' & '
    equation_string += '='
    equation_string += bmatrix(rhs,'b')
    equation_string += '\end{matrix}$'
    return equation_string

###############################################################################
# STREAMLIT
###############################################################################

st.set_page_config(layout="wide")

# setup session_states
if 'selected_member' not in st.session_state:
    st.session_state.selected_member = []
if 'selected_nodes' not in st.session_state:
    st.session_state.selected_nodes = []

st.title("Calculating internal forces of a beam structure")

###############################################################################
# INPUTS
###############################################################################

st.sidebar.write('you need to fulfill 2*n_nodes = n_members + n_supports to get a square matrix')

nodes_str = st.sidebar.text_input(label = "nodes", help = "[x-position,y-position]", value='''[0,0],
[1,1],
[1,0],
[2,2],
[2,0],
[3,1],
[3,0],
[4,0]''')

exec("nodes = np.array([" + nodes_str + "])")

members_str = st.sidebar.text_input(label = "members", help = "[1st node,2nd node]", value='''[1,2],
[1,3],
[2,3],
[3,5],
[2,5],
[2,4],
[4,5],
[5,7],
[5,6],
[6,7],
[4,6],
[7,8],
[6,8]''')

exec("members = np.array([" + members_str + "])")

if not 'members' in st.session_state:
    st.session_state.members = np.array([[1,2],
        [1,3],
        [2,3],
        [3,5],
        [2,5],
        [2,4],
        [4,5],
        [5,7],
        [5,6],
        [6,7],
        [4,6],
        [7,8],
        [6,8]
        ])
members = st.session_state.members

support_str = st.sidebar.text_input(label = "support", help = "[node,angle]", value='''[1, 0],
[1, 90],
[8, 90]''')

exec("support = np.array([" + support_str + "],dtype='f')")

f_ext_str = st.sidebar.text_input(label = "external forces", help = "[node,angle,force]", value='''[3,-90,10],
[4,180,10],
[5,-90,15]''')

exec("f_ext = np.array([" + f_ext_str + "],dtype='f')")

# convert angles to radians
support[:,1] = np.radians(support[:,1])
f_ext[:,1] = np.radians(f_ext[:,1])

###############################################################################
# SIDEBARS
###############################################################################

decimals = st.sidebar.number_input(label="precision of print",min_value=0,max_value=5,value=2)
textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize',r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)

###############################################################################
# CALCULATIONS
###############################################################################

matrix, rhs, internal_forces = update_data(nodes,members,support,f_ext)

[fig,forcemap] = update_plot(internal_forces,members,nodes,f_ext,support)

#                 np.zeros([2*len(nodes),1])
# for i in np.arange(0,len(label_vecotr)):
#     label_vector

###############################################################################
# OUTPUTS
###############################################################################

with st.expander('Look at the plot', expanded=True):
    sn = plotly_events(fig)#, click_event=True)
    #st.plotly_chart(forcemap)
    st.sidebar.write('return value of plotly_events: ' + str(sn))
    if not sn == []:
        if sn[0]['curveNumber'] == len(members)+0:
            st.session_state.selected_nodes.append(sn[0]['pointNumber'])
            st.write('You selected node #'
                         + str(sn[0]['pointNumber']+1)
                         + '. Select another one to draw a new beam')
            st.write('current st.session_state.selected_nodes (actual python inidces, mind you, not the fancy n+1 indices): ' + str(st.session_state.selected_nodes))
            if len(st.session_state.selected_nodes) == 2:
                #st.session_state.selected_nodes[1] = sn[0]['pointNumber']
                #new_member(st.session_state.selected_nodes)
                st.session_state.selected_nodes = []
            elif len(st.session_state.selected_nodes) > 2:
                st.session_state.selected_nodes = []
            
                
                
        if sn[0]['curveNumber'] == len(members)+1:
            st.write('You selected a force. Forces can only be set in the sidebar.')
        if sn[0]['curveNumber'] == len(members)+2:
            st.write('You selected a support. Supports can only be set in the sidebar.')
        if sn[0]['curveNumber'] == len(members)+3:
            st.session_state.selected_member.append(sn[0]['pointNumber'])

with st.expander('Look at the Matrix. Select font size in the sidebar', expanded=True):
    st.markdown(print_equations(matrix, rhs, internal_forces,len(members),len(support),decimals,textsize))