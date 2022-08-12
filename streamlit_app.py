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

#import math

###############################################################################
# TO-DO
###############################################################################


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
    #st.write(str(connected_nodes))
    forces_connected = True
    for f in force_range:
        if int(f_ext[f,0]) in deletenodes:
            forces_connected = False
    
    issquare = 2*len(connected_nodes) == (len(state.members) + len(state.support))
    
    return rods_per_node, connected_nodes, issquare, forces_connected

# maybe use @st.cache()?
def update_data(nodes,members,support,f_ext,debug,rods_per_node,buildonly=False):
    
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
    matrix = np.append(matrix, np.zeros([len(matrix),len(support)]),axis = 1)
    for s in support_range:
        inode = int(support[s,0])
        angle = np.radians(support[s,1])
        if debug:
            st.write('writing support #' + str(s) + ' at node #' + str(inode) + ' to matrix')
        matrix[2*inode,n_members+s] = np.cos(angle) # forces along x
        matrix[2*inode+1,n_members+s] = np.sin(angle) # forces along y
        
    # compute right hand side
    rhs = np.zeros([2*n_nodes,1])
    for force in force_range:
        inode = int(f_ext[force,0])
        angle = np.radians(f_ext[force,1])
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
    if debug:
        #connected_nodes=np.delete(nodes,deletenodes,0)
        st.write('nodes: ' + str(nodes) + ' n_nodes: ' + str(n_nodes))
        st.markdown(print_equations(matrix, rhs, [],len(members),len(state.support),1,'\scriptsize'))
        st.write('connected nodes: ' + str(connected_nodes))
        st.write('delete rows: ' + str(deleterows))
    matrix = np.delete(matrix,deleterows,0)
    rhs = np.delete(rhs,deleterows,0)    
    #st.write(str(matrix))
        
    if debug:
        st.markdown(print_equations(matrix, rhs, [],len(members),len(state.support),1,'\scriptsize'))
        
    # solve for unknowns
    if buildonly:
        internal_forces = []
    else:
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

def new_force_input():
    if not state.new_f_ext_str == '':
        new_f_ext_list = string_to_list(state.new_f_ext_str)
        #new_f_ext_list[1] = math.radians(new_f_ext_list[1])
        if state.new_f_ext == [[]]:
            state.new_f_ext = [new_f_ext_list]
        else:
            np.append(state.new_f_ext,[new_f_ext_list],axis=0)
        state.new_f_ext_str = ''
    return

def new_support_input():
    if not state.new_support_str == '':
        new_support_list = string_to_list(state.new_support_str)
        #new_support_list[1] = math.radians(new_support_list[1])
        if state.new_supports == [[]]:
            state.new_supports = [new_support_list]
        else:
            if new_support_list not in state.new_supports:
                state.new_supports = np.append(state.new_supports,[new_support_list],axis=0)
        state.new_support_str = ''
    return

def new_force_delete():
    if not state.force_delete_str == '':
        if len(state.force_delete_str) > 2:
            st.warning('Delete one force at a time!')
        else:
            iforce = int(state.force_delete_str)
            if iforce not in state.removed_forces:
                state.removed_forces.append(iforce)
        state.force_delete_str = ''
    return

def new_support_delete():
    if not state.support_delete_str == '':
        if len(state.support_delete_str) > 2:
            st.warning('Delete one support at a time!')
        else:
            isupport = int(state.support_delete_str)
            if isupport not in state.removed_supports:
                state.removed_supports.append(isupport)
        state.support_delete_str = ''
    return

def new_array(arraytype):
    if arraytype == 'support':
        exec("state.support = np.array([" + state.support_str + "])")
    elif arraytype == 'f_ext':
        exec("state.f_ext = np.array([" + state.f_ext_str + "])")
    elif arraytype == 'members':
        exec("state.members = np.array([" + state.members_str + "])")
    return

def string_to_list(stringlist):
    list_of_str = stringlist.split()
    list_from_str = [float(x) for x in list_of_str]
    return list_from_str

def update_all(removed_members,new_members,removed_forces,removed_supports,new_f_ext,new_supports):
    if not removed_members == []:
        state.members = np.delete(state.members,removed_members,axis=0)
        state.removed_members = []
    if not new_members == [[]]:
        state.members = np.append(state.members,new_members,axis=0)
        state.new_members = [[]]
        
    if not removed_forces == []:
        state.f_ext = np.delete(state.f_ext,removed_forces,axis=0)
        state.removed_forces = []
    if not new_f_ext == [[]]:
        if debug:
            st.write(str(new_f_ext))
        state.f_ext = np.append(state.f_ext,new_f_ext,axis=0)
        state.new_f_ext = [[]]
        
    if not removed_supports == []:
        state.support = np.delete(state.support,removed_supports,axis=0)
        state.removed_supports = []
    if not new_supports == [[]]:
        state.support = np.append(state.support,new_supports,axis=0)
        state.new_supports = [[]]
        
    return

def reset_data():
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
    state.support = np.array([
        [0, 0],[0, 90],[20, 90]
        ])
    state.f_ext = np.array([
        [5,-90,10],[12,180,10],[15,-90,15]
        ])
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

def update_plot(internal_forces,members,nodes,f_ext,support,onlyviz):

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
        if onlyviz:
            fig.add_trace(go.Scatter(
                x=x,
                y=y,
                mode="lines + markers",
                name='member #' + str(m),
                showlegend=False,
                line=dict(
                    color='black'
                    )
                ))
        else:
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
        angle = np.radians(f_ext[f,1])
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
        angle = np.radians(support[s,1])
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
        if onlyviz:
            member_names.append('deselect member #' + str(m))
        else:
            member_names.append('member #' + str(m) + '; current force: ' + str(round(internal_forces[m][0])))
    if onlyviz:    
        fig.add_trace(go.Scatter(
            x=centresx,
            y=centresy,
            mode='markers',
            marker=dict(symbol='x',color='orange'),
            name='members for deselection',
            text=member_names
            ))
    else:
        fig.add_trace(go.Scatter(
            x=centresx,
            y=centresy,
            mode='markers',
            marker=dict(symbol='circle',color='orange'),
            name='',
            text=member_names
            ))
        
    return fig,forcemap

###############################################################################
# PRINTING OUTPUT
###############################################################################

def bmatrix(a,showzeros,matrixtype=''):
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
            if a[x][y] == 0:
                text += showzeros
            else:
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

def latex_to_md(textsize):
    string = r'''<font size="'''
    if textsize == r'\tiny':
        string+='1'
    elif textsize == r'\scriptsize':
        string+='2'
    elif (textsize == r'\footnotesize') | (textsize == r'\small'):
        string+='3'
    elif textsize == r'\normalsize':
        string+='4'
    string += r'''"> '''
    return string

def print_equations(matrix, rhs, internal_forces,n_rods,n_bcs,decimals,textsize,n_nodes,onlyviz,showzeros):
    label_vector = []
    k=0
    for node in connected_nodes:
        label_vector.append([r' \text{'
                             #+ 'N.' + str(int(np.floor(k/2))) + ' '
                             + str(node)])
        label_vector[k][0] += ', x'
        label_vector[k][0] += '}'
        k+=1
        label_vector.append([r' \text{'
                             #+ 'N.' + str(int(np.floor(k/2))) + r' '
                             + str(node)])
        label_vector[k][0] += ', y'
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
    equation_string += r'\text{' + str(n_nodes) + r' nodes} & \text{' + str(n_rods) + ' rods and ' + str(n_bcs) + ' boundary conditions}'
    #equation_string += bmatrix(label_vector,'h')
    equation_string += r' & \text{internal forces} & \text{external forces}\\'
    equation_string += bmatrix(label_vector,showzeros,'v')
    equation_string += ' & '
    equation_string += bmatrix(matrix,showzeros,'b')
    equation_string += ' & '
    equation_string += '\cdot '
    if (onlyviz) | (internal_forces == []):
        equation_string += '?'
    else:
        equation_string += bmatrix(internal_forces,showzeros,'b')
    equation_string += ' & '
    equation_string += '='
    equation_string += bmatrix(rhs,showzeros,'b')
    equation_string += '\end{matrix}$'
    return equation_string

###############################################################################
# STREAMLIT AND STATES
###############################################################################

st.set_page_config(layout="wide",initial_sidebar_state='collapsed')

# setup session_states
if 'removed_members' not in state:
    state.removed_members = []
if 'new_members' not in state:
    state.new_members = [[]]
if 'selected_nodes' not in state:
    state.selected_nodes = []
if 'members_str' not in state:
    state.members_str = ''
    
if 'removed_forces' not in state:
    state.removed_forces = []
if 'new_f_ext' not in state:
    state.new_f_ext = [[]]
if 'new_f_ext_str' not in state:
    state.new_f_ext_str = ''
if 'f_ext_str' not in state:
    state.f_ext_str = ''
if 'force_delete_str' not in state:
    state.force_delete_str = ''
    
if 'removed_supports' not in state:
    state.removed_supports = []
if 'new_supports' not in state:
    state.new_supports = [[]]
if 'new_support_str' not in state:
    state.new_support_str = ''
if 'support_str' not in state:
    state.support_str = ''
if 'support_delete_str' not in state:
    state.support_delete_str = ''

# if 'internal_forces' not in state:
#     state.internal_forces = np.array([])
# if 'matrix' not in state:
#     state.matrix = []
# if 'rhs' not in state:
#     state.rhs = []
    
if 'all_nodes' not in state:
    l = 5#st.sidebar.number_input(label='length of plot',min_value=2,max_value=10,value=5)
    h = 4#st.sidebar.number_input(label='height of plot',min_value=1,max_value=10,value=4)
    all_nodes = []
    for x in range(l+1):
        for y in range(h+1):
            all_nodes.append([x,y])
    state.all_nodes=np.array(all_nodes)

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
    
st.title("Internal forces of a rod structure")

with st.expander('Explanation'):
    st.write("Here, you can see a rod system with supports and forces. Below, you see the linear system of equations that is solved to calculate the iternal forces of the rods.  \n "
             + "  \n "
             + "Select 'Interactive mode' to be able to change the rod structure.  \n "
             + "Delete rods via the orange 'x' markers. Add rods by klicking both nodes after each other (you may need to hide the rod-markers by klicking on 'x members for deselection' in the legend). Below the plot, you get information about which node(s) you chose.  \n "
             + "Add or delete supports or forces in the sidebar.  \n "
             + "You can input all forces, supports and members in bulk by choosing 'Use vectors as input' in the sidebar.  \n "
             + "  \n "
             + "To display your changes, klick 'Update plot and matrix'.  \n "
             + "To update your calculations, click 'Update calculations' and deselect 'Interactive mode'.  \n "
             + "  \n "
             + "Keep in mind that your rod structure must consist of triangles to yield sensible results!  \n "
             )

###############################################################################
# DEBUG OPTIONS
###############################################################################

debug = False#st.sidebar.checkbox(label="show development stuff")

if debug:
    decimals = st.sidebar.number_input(label="precision of print",min_value=0,max_value=5,value=2)
    #textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize',r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)
else:
    decimals = 2
    #textsize = r'\scriptsize'

if not debug:
    if 'support' not in state:
        state.support = np.array([
            [0, 0],[0, 90],[20, 90]
            ])
    if 'f_ext' not in state:
        state.f_ext = np.array([
            [5,-90,10],[12,180,10],[15,-90,15]
            ])

###############################################################################
# BUTTONS
###############################################################################

col1,col2,col3,col4 = st.columns(4)#([1,1,1,1,0])
with col1:
    reset = st.button('Reset data')
if reset: reset_data()
with col2:
    onlyviz = st.checkbox("Interactive mode",key='onlyviz')
with col3:
    apply_changes = st.button('Update plot and matrix')#,on_click=update_members(state.removed_members,state.new_members))
with col4:
    calculate = st.button('Update calculations')
    if calculate: onlyviz = False

# if onlyviz:
#     st.write('Mind the admonitions below the plot about requirements to the rod system.')

###############################################################################
# SIDEBARS
###############################################################################

textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize',r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)
textsize_md = latex_to_md(textsize)

showzeros = st.sidebar.selectbox(label="show zeros as...", options=[' 0 ',' '])

if onlyviz:
    vectorinput = st.sidebar.checkbox("Use vectors as input. (klick 'Update plot' once to initialize)")
    if vectorinput:
        members_str = st.sidebar.text_input(label = "members", 
                                            help = "[node,angle]", 
                                            value='''[0,6],[0,5],[6,5],[5,10],[6,10],[6,12],[12,10],[10,15],[10,16],[16,15],[12,16],[15,20],[16,20]''',
                                            key='members_str',
                                            on_change=new_array('members'))
        
        support_str = st.sidebar.text_input(label = "support", 
                                            help = "[node,angle]", 
                                            value='[0, 0],[0, 90],[20, 90]',
                                            key='support_str',
                                            on_change=new_array('support'))
        
        f_ext_str = st.sidebar.text_input(label = "external forces", 
                                          help = "[node,angle,force]", 
                                          value='[5,-90,10],[12,180,10],[15,-90,15]',
                                          key='f_ext_str',
                                          on_change=new_array('f_ext'))
        
    else:
        new_f_ext_str = st.sidebar.text_input(label = "add external forces", 
                                          help = "node angle force",
                                          placeholder="e.g. '6 90 10'", 
                                          on_change=new_force_input(), 
                                          key='new_f_ext_str')
        
        new_support_str = st.sidebar.text_input(label = "add supports", 
                                            help = "node angle",
                                            placeholder="e.g. '10 90'", 
                                            on_change=new_support_input(),
                                            key='new_support_str')
        
        f_ext_text = ''
        for i in range(len(state.f_ext)):
                f_ext_text += str(i) + ': ' + str(state.f_ext[i])
                if i < len(state.f_ext)-1:
                        f_ext_text += ', '
        support_text = ''
        for i in range(len(state.support)):
                support_text += str(i) + ': ' + str(state.support[i])
                if i < len(state.support)-1:
                        support_text += ', '
        
        st.sidebar.markdown(textsize_md + 'Currently, you have these forces: '
                            + f_ext_text + '. Would you like to delete one? '
                            + r''' </font>''', unsafe_allow_html=True)
        st.sidebar.text_input(label='remove external forces', 
                              placeholder="IDs start at 0!",
                              on_change=new_force_delete(),
                              key='force_delete_str')
        
        st.sidebar.markdown(textsize_md + 'Currently, you have these supports: '
                            + support_text + '. Would you like to delete one? '
                            + r''' </font>''', unsafe_allow_html=True)
        st.sidebar.text_input(label='remove supports',
                              placeholder="IDs start at 0!",
                              on_change=new_support_delete(),
                              key='support_delete_str')
        
        if apply_changes:
            update_all(state.removed_members,state.new_members,state.removed_forces,state.removed_supports,state.new_f_ext,state.new_supports)
            state.selected_nodes = []
          
###############################################################################
# CALCULATIONS
###############################################################################

[rods_per_node, connected_nodes, issquare, forces_connected] = check_data(state.all_nodes,state.members,state.support,state.f_ext)

if 'matrix' not in state:
    [state.matrix, state.rhs, state.internal_forces] = update_data(state.all_nodes,state.members,state.support,state.f_ext,debug,rods_per_node)
if 'fig' not in state:
    [state.fig,state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)

if onlyviz:
    if state.onlyviz:
        [state.fig,state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
    if apply_changes:
        [state.fig,state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
elif calculate:
    if issquare:
        [state.matrix, state.rhs, state.internal_forces] = update_data(state.all_nodes,state.members,state.support,state.f_ext,debug,rods_per_node)
        [state.fig,state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
    else:
        st.warning("I am having issues solving your system. Return to interactive mode to check whether your system is solvable. If you need to reset, press 'Reste data'.")
        onlyviz = True
        [state.fig,state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
    
###############################################################################
# OUTPUTS
###############################################################################

# loading fig from cache
fig = state.fig

if onlyviz:
    sn = plotly_events(fig)
    if debug:
        st.write('return value of plotly_events: ' + str(sn))
    if not sn == []:
        
        if sn[0]['curveNumber'] == len(state.members)+0:
            state.selected_nodes.append(sn[0]['pointNumber'])
            if len(state.selected_nodes) == 1:
                st.markdown(textsize_md + 'You selected node #'
                         + str(sn[0]['pointNumber'])
                         + '. Select another one to draw a new rod. '
                         + r''' </font>''', unsafe_allow_html=True)
            elif len(state.selected_nodes) == 2:
                if not state.selected_nodes[1] == state.selected_nodes[0]:
                    new_member(state.selected_nodes)
                    st.markdown(textsize_md + 'You selected nodes #'
                         + str(state.new_members[-1])
                         + ". See below which changes will be applied by 'update plot'. "
                         + r''' </font>''', unsafe_allow_html=True)
                    state.selected_nodes = []
                else:
                    st.markdown(textsize_md + 'You selected node #'
                             + str(sn[0]['pointNumber'])
                             + '. Twice. Select a different one. '
                             + r''' </font>''', unsafe_allow_html=True)
            elif len(state.selected_nodes) > 2:
                state.selected_nodes = []
                st.markdown(textsize_md + 'You selected more than 2 nodes. Must be a bug. Cleared selected nodes.'
                         + r''' </font>''', unsafe_allow_html=True)
            if debug:
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
                
        # Will deselect supports and forces explicitly until I find out how exactly the nodes in the quiver are numbered    
        # if sn[0]['curveNumber'] == len(state.members)+1:
        #     if (sn[0]['pointNumber']-10)%4 == 0:
        #         iforce = int((sn[0]['pointNumber']-10)/4)
        #         if iforce not in state.removed_forces:
        #             state.removed_forces.append(iforce)
        #     else:
        #         st.write('Please select the force via the tip of the arrow.')
            
        # if sn[0]['curveNumber'] == len(state.members)+2:
        #     if (sn[0]['pointNumber']-10)%4 == 0:
        #         isupport = int((sn[0]['pointNumber']-10)/4)
        #         if isupport not in state.removed_supports:
        #             state.removed_supports.append(isupport)
        #     else:
        #         st.write('Please select the support via the tip of the arrow.')
            
        if sn[0]['curveNumber'] == len(state.members)+3:
            if sn[0]['pointNumber'] not in state.removed_members:
                state.removed_members.append(sn[0]['pointNumber'])
            
    # checking if we will get a square matrix
    if forces_connected == False:
        st.warning('Not all forces are connected. You may solve the system anyway.')
    
    st.markdown(textsize_md
                + r"'Update plot' will remove **members** #" + str(state.removed_members)
                + ' and add members on ' + str(state.new_members) + '. '
                + r'You remove **supports** #' + str(state.removed_supports)
                + 'and add supports on ' + str(np.round(state.new_supports,2)) + '. '
                + r'You remove **forces** #' + str(state.removed_forces)
                + ' and add forces on ' + str(np.round(state.new_f_ext,2)) + '. '
                + r''' </font>''', unsafe_allow_html=True)
    
    status_string = (textsize_md + r'''You need to fulfill 
        $$ 2 \cdot n_\text{nodes} = n_\text{members} + n_\text{supports} $$
        to get a square matrix. ''' + 'You have '
        + str(len(connected_nodes)) + ' nodes, '
        + str(len(state.members)) + ' members and '
        + str(len(state.support)) + ' supports. ')
    if issquare:
        status_string += 'Currently, you get a square matrix.'
    else:
        if (2*len(connected_nodes) > (len(state.members) + len(state.support))):
            status_string += 'You need to add more members or supports.'
        else:
            status_string += 'You need to remove members or supports or add more nodes.'
    st.markdown(status_string + r''' </font>''', unsafe_allow_html=True)
        
else:
    st.plotly_chart(fig)

if onlyviz:
    buildonly = True
    matrix, rhs, internal_forces = update_data(state.all_nodes,state.members,state.support,state.f_ext,debug,rods_per_node,buildonly)
    st.markdown(
        print_equations(
            matrix, rhs, internal_forces,
            len(state.members),len(state.support),
            decimals,textsize,len(connected_nodes),
            onlyviz,showzeros))
else:
    st.markdown(
        print_equations(
            state.matrix, state.rhs, state.internal_forces,
            len(state.members),len(state.support),
            decimals,textsize,len(connected_nodes),
            onlyviz,showzeros))
    
