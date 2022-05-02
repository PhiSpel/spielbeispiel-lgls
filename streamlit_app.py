# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:38:24 2022

@author: Philipp Spelten
"""

import numpy as np
# import pandas as pd

# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator
# import matplotlib.colors as mcol
# import matplotlib.cm as cm

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
from streamlit import session_state as state
from streamlit_plotly_events import plotly_events

# import math

###############################################################################
# TO-DO
###############################################################################

# 3. create error-messages

###############################################################################
# CALCULATIONS
###############################################################################

def check_data(nodes,members,support,f_ext):
    n_members = len(members)
    members_range = np.arange(0,n_members)
    n_nodes = len(nodes)
    nodes_range = np.arange(0,n_nodes)
    nsupport = len(support)
    support_range = np.arange(0,nsupport)
    nf_ext = len(f_ext)
    force_range = np.arange(0,nf_ext)
    rods_per_node = np.zeros(len(nodes))
    for m in members_range:
        node1 = members[m,0]
        node2 = members[m,1]
        rods_per_node[node1] += 1
        rods_per_node[node2] += 1
    for s in support_range:
        inode = int(support[s,0])
        rods_per_node[inode] += 1
    for f in force_range:
        inode = int(f_ext[f,0])
        rods_per_node[inode] += 1
    deletenodes = []
    for inode in nodes_range:
        if rods_per_node[inode] < 2:
            deletenodes.append(inode)
    connected_nodes=np.delete(nodes,deletenodes,0)
    for f in force_range:
        if int(f_ext[f,0]) not in connected_nodes:
            forces_connected = False
    
    issquare = 2*len(connected_nodes) == (len(state.members) + len(state.support))
    
    return rods_per_node, connected_nodes, issquare, forces_connected

# maybe use @st.cache()?
def update_data(nodes,members,support,f_ext,debug,rods_per_node):
    
    n_members = len(members)
    members_range = np.arange(0,n_members)
    n_nodes = len(nodes)
    nodes_range = np.arange(0,n_nodes)
    nsupport = len(support)
    support_range = np.arange(0,nsupport)
    nf_ext = len(f_ext)
    force_range = np.arange(0,nf_ext)
    
    members = np.append(members, np.zeros([n_members,3]),axis = 1)
    # compute angles
    for m in members_range:
        node1 = int(members[m,0])
        node2 = int(members[m,1])
        stab = nodes[node2,:]-nodes[node1,:]
        direction = np.sign(stab[1])
        if direction == 0:
            direction = 1
        alpha = np.arccos(stab[0]/np.linalg.norm(stab))*direction
        members[m,2] = np.cos(alpha) # along x
        members[m,3] = np.sin(alpha) # along y
        members[m,4] = alpha
    
    # compute matrix coefficients
    ks1 = np.zeros([n_nodes,n_members])
    ks2 = np.zeros([n_nodes,n_members])
    matrix = np.zeros([2*n_nodes,n_members])
    for k in nodes_range:
        # compute x and y forces per member
        for m in members_range:
            node1 = int(members[m,0])
            node2 = int(members[m,1])
            ks1[k,m] = int(k==node1)
            ks2[k,m] = -int(k==node2)
        ks = ks1+ks2
        # forces along x
        matrix[2*k,:] = np.multiply(ks[k,:],members[:,2])
        # forces along y
        matrix[2*k+1,:] = np.multiply(ks[k,:],members[:,3])
        
    # add support forces to matrix
    if debug:
        st.write('supports: ' + str(support))
    matrix = np.append(matrix, np.zeros([len(matrix),3]),axis = 1)
    for s in support_range:
        inode = int(support[s,0])
        angle = support[s,1]
        matrix[2*inode,n_members+s] = np.cos(angle) # forces along x
        matrix[2*inode+1,n_members+s] = np.sin(angle) # forces along y
        
    # compute right hand side
    rhs = np.zeros([2*n_nodes,1])
    for force in force_range:
        inode = int(f_ext[force,0])
        angle = f_ext[force,1]
        newtons = f_ext[force,2]
        # external force along x
        rhs[2*inode,0] = newtons*np.cos(angle)
        # external force along y
        rhs[2*inode+1,0] = newtons*np.sin(angle)
        
    # delete obsolete rows
    deleterows = []
    for inode in nodes_range: #number of the node we are talking about
        k = 2*inode
        if rods_per_node[inode] <= 1:
            deleterows.append(k)
            deleterows.append(k+1)
    # if debug:
    #     connected_nodes=np.delete(nodes,deletenodes,0)
    #     st.write('nodes: ' + str(nodes) + ' n_nodes: ' + str(n_nodes))
    #     st.markdown(print_equations(matrix, rhs, [],len(members),len(state.support),1,'\scriptsize'))
    #     st.write('connected nodes: ' + str(connected_nodes))
    #     st.write('delete rows: ' + str(deleterows))
    matrix = np.delete(matrix,deleterows,0)
    rhs = np.delete(rhs,deleterows,0)    
    #st.write(str(matrix))
        
    if debug:
        st.markdown(print_equations(matrix, rhs, [],len(members),len(state.support),1,'\scriptsize'))
        
    # solve for unknowns
    internal_forces = np.linalg.solve(matrix, rhs)
    
    return matrix, rhs, internal_forces

def new_member(new_nodes):
    node1 = new_nodes[0]
    node2 = new_nodes[1]
    if state.new_members == [[]]:
        state.new_members = [[node1,node2]]
    else:
        state.new_members = np.append(state.new_members,[[node1,node2]],axis=0)
    return

def update_members(removed_members,new_members):
    if not ((new_members == [[]]) & (removed_members == [])):
        state.members = np.delete(state.members,removed_members,axis=0)
        state.members = np.append(state.members,new_members,axis=0)
        state.new_members = [[]]
        state.removed_members = []
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
    
    # draw rods
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
        node1_index = int(members[m,0])
        node2_index = int(members[m,1])
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
            name='member #' + str(m),
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
                    text=np.arange(0,len(nodes)),
                    name='nodes'
                    ))
    
    # calculate coordinates of forces
    x0 = []
    y0 = []
    fx = []
    fy = []
    for f in np.arange(0,len(f_ext)):
        node_index = int(f_ext[f,0])
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
    support_length = 5
    x0 = []
    y0 = []
    sx = []
    sy = []
    for s in np.arange(0,len(support)):
        node_index = int(support[s,0])
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
        
    #draw center points to be able to deselect rods
    centresx=[]
    centresy=[]
    member_names=[]
    for m in np.arange(0,len(members)):
        node1_index = int(members[m,0])
        node2_index = int(members[m,1])
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

def print_equations(matrix, rhs, internal_forces,n_rods,n_bcs,decimals,textsize):
    label_vector = []
    k=0
    for node in connected_nodes:
        label_vector.append([r' \text{' + str(node)])
        label_vector[k][0] += 'x'
        label_vector[k][0] += '}'
        k+=1
        label_vector.append([r' \text{' + str(node)])
        label_vector[k][0] += 'y'
        label_vector[k][0] += '}'
        k+=1
    label_vector=np.array(label_vector)
    matrix = np.round(matrix,decimals)
    matrix[matrix==0] = 0
    rhs = np.round(rhs,decimals)
    rhs[rhs==0] = 0
    internal_forces = np.round(internal_forces,decimals)
    internal_forces[internal_forces==0] = 0
    equation_string = r'$'
    equation_string += textsize
    equation_string += r' \begin{matrix} '
    equation_string += r'\text{nodes} & \text{' + str(n_rods) + ' rods and ' + str(n_bcs) + ' boundary conditions}'
    #equation_string += bmatrix(label_vector,'h')
    equation_string += r' & \text{internal forces} & \text{external forces}\\'
    equation_string += bmatrix(label_vector,'v')
    equation_string += ' & '
    equation_string += bmatrix(matrix,'b')
    equation_string += ' & '
    equation_string += '\cdot '
    if not internal_forces == []:
        equation_string += bmatrix(internal_forces,'b')
    else:
        equation_string += '?'#bmatrix(internal_forces,'b')
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
if 'removed_members' not in state:
    state.removed_members = []
if 'new_members' not in state:
    state.new_members = [[]]
if 'selected_nodes' not in state:
    state.selected_nodes = []

st.title("Internal forces of a rod structure")

###############################################################################
# INPUTS
###############################################################################

l = 5#st.sidebar.number_input(label='length of plot',min_value=2,max_value=10,value=5)
h = 4#st.sidebar.number_input(label='height of plot',min_value=1,max_value=10,value=4)

all_nodes = []
for x in range(l+1):
    for y in range(h+1):
        all_nodes.append([x,y])
all_nodes=np.array(all_nodes)

if 'members' not in state:
    state.members = np.array([
        [0,6],
        [0,5],
        [6,5],
        [5,10],
        [6,10],
        [6,12],
        [12,10],
        [10,15],
        [10,16],
        [16,15],
        [12,16],
        [15,20],
        [16,20]
        ])
    
# if 'support' not in state:
#     state.support = np.array([
#         [0, 0],
#         [0, 90],
#         [20, 90]
#         ],dtype='f')
#     # convert angles to radians
#     state.support[:,1] = np.radians(state.support[:,1])
    
# if 'f_ext' not in state:
#     state.f_ext = np.array([
#         [5,-90,10],
#         [12,180,10],
#         [10,-90,15]
#         ],dtype='f')
#     # convert angles to radians
#     state.f_ext[:,1] = np.radians(state.f_ext[:,1])

support_str = st.sidebar.text_input(label = "support", help = "[node,angle]", 
                                    value='''[0, 0],[0, 90],[20, 90]''')

if 'support' not in state: 
    exec("state.support = np.array([" + support_str + "],dtype='f')")
    state.support[:,1] = np.radians(state.support[:,1])
# else:
#     exec("state.support = np.array([" + support_str + "],dtype='f')")

f_ext_str = st.sidebar.text_input(label = "external forces", help = "[node,angle,force]", value='''[5,-90,10],[12,180,10],[15,-90,15]''')

if 'f_ext' not in state: 
    exec("state.f_ext = np.array([" + f_ext_str + "],dtype='f')")
    state.f_ext[:,1] = np.radians(state.f_ext[:,1])
# else:
#     exec("state.f_ext = np.array([" + f_ext_str + "],dtype='f')")
   
st.sidebar.write('members: ' + str(state.members))
#st.sidebar.write('nodes: ' + str(all_nodes))
#st.sidebar.write('support: ' + str(state.support.round()))
#st.sidebar.write('f_ext: ' + str(state.f_ext.round()))

if 'internal_forces' not in state:
    state.internal_forces = []
if 'matrix' not in state:
    state.matrix = []
if 'rhs' not in state:
    state.rhs = []

###############################################################################
# SIDEBARS
###############################################################################

debug = st.sidebar.checkbox(label="show development stuff")
if debug:
    decimals = st.sidebar.number_input(label="precision of print",min_value=0,max_value=5,value=2)
    textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize',r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)
else:
    decimals = 2
    textsize = r'\scriptsize'
###############################################################################
# VISUAL VS CALCULATED
###############################################################################

[col1,col2] = st.columns([4,1])

with col1:
    onlyviz = st.checkbox("Visualization only. Choose this to be able to change the structure.")

with col2:
    if onlyviz:
        apply_changes = st.button('Apply changes')#,on_click=update_members(state.removed_members,state.new_members))
        if apply_changes:
            update_members(state.removed_members,state.new_members)
            
if onlyviz:
    st.write('Mind the admonitions below the plot on requirements to the rod system.')

###############################################################################
# CALCULATIONS
###############################################################################

[rods_per_node, connected_nodes, issquare, forces_connected] = check_data(all_nodes,state.members,state.support,state.f_ext)

if onlyviz:
    # calculating
    [fig,forcemap] = update_plot(state.internal_forces,state.members,all_nodes,state.f_ext,state.support)
else:
    if issquare:
        state.matrix, state.rhs, state.internal_forces = update_data(all_nodes,state.members,state.support,state.f_ext,debug,rods_per_node)
        [fig,forcemap] = update_plot(state.internal_forces,state.members,all_nodes,state.f_ext,state.support)
    else:
        st.warning("I am having issues solving your system. Return to visualization only to check whether your system is solvable. If you neet to reset, press 'c' or refresh your browser.")
        [fig,forcemap] = update_plot(state.internal_forces,state.members,all_nodes,state.f_ext,state.support)
    

###############################################################################
# OUTPUTS
###############################################################################

if onlyviz:
    with st.expander('Look at the plot', expanded=True):
        sn = plotly_events(fig)
        if debug:
            st.write('return value of plotly_events: ' + str(sn))
        if not sn == []:
            if sn[0]['curveNumber'] == len(state.members)+0:
                state.selected_nodes.append(sn[0]['pointNumber'])
                if len(state.selected_nodes) == 1:
                    st.write('You selected node #'
                             + str(sn[0]['pointNumber'])
                             + '. Select another one to draw a new rod')
                nodes_string = 'Current storage in state.selected_nodes: ' + str(state.selected_nodes)
                if len(state.selected_nodes) == 2:
                    if state.selected_nodes[1] == state.selected_nodes[0]:
                        nodes_string += '. Will not write to state.new_members, b/c you selected the same node twice.'
                    else:
                        new_member(state.selected_nodes)
                        nodes_string += '. Wrote nodes to state.new_members.'
                    state.selected_nodes = []
                elif len(state.selected_nodes) > 2:
                    state.selected_nodes = []
                    nodes_string += '. Must be a bug. Cleared selected nodes.'
                st.write(nodes_string)
                
            if sn[0]['curveNumber'] == len(state.members)+1:
                st.write('You selected a force. Forces can only be set in the sidebar.')
            if sn[0]['curveNumber'] == len(state.members)+2:
                st.write('You selected a support. Supports can only be set in the sidebar.')
            if sn[0]['curveNumber'] == len(state.members)+3:
                state.removed_members.append(sn[0]['pointNumber'])
    # checking if we will get a square matrix
    if forces_connected == False:
        st.warning('Not all forces are connected. You may solve the system anyway.')
    status_string = r'''You need to fulfill $$ 2 \cdot n_\text{nodes} = n_\text{members} + n_\text{supports} $$ to get a square matrix. '''
    status_string += 'You have ' + str(len(connected_nodes)) + ' nodes, ' + str(len(state.members)) + ' members and ' + str(len(state.support)) + ' supports. '
    if issquare:
        status_string += 'Currently, you get a square matrix.'
    else:
        if (2*len(connected_nodes) > (len(state.members) + len(state.support))):
            status_string += 'You need to add more members or supports.'
        else:
            status_string += 'You need to add more nodes.'
    st.write(status_string)
    st.markdown(r'You chose to remove members #' + str(state.removed_members) + r''' You chose to add members on ''' + str(state.new_members))
else:
    st.plotly_chart(fig)

with st.expander('Look at the Matrix. Select font size in the sidebar', expanded=True):
    st.markdown(print_equations(state.matrix, state.rhs, state.internal_forces,len(state.members),len(state.support),decimals,textsize))
    
