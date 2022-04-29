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
from plotly.figure_factory import create_quiver

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

def update_plot(internal_forces,members,nodes,f_ext):

    fig = go.Figure()
    
    # draw nodes
    fig.add_trace(go.Scatter(x=nodes[:,0],y=nodes[:,1],
                    mode='markers',
                    text=np.arange(1,len(nodes)+1)))
    
    # Make a user-defined colormap.
    #cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","r","w","b","b"])
    # Make a normalizer that will map the time values from
    # [start_time,end_time+1] -> [0,1].
    #cnorm = mcol.Normalize(vmin=min(internal_forces),vmax=max(internal_forces))
    # Turn these into an object that can be used to map time values to colors and
    # can be passed to plt.colorbar().
    #cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    #cpick.set_array([])
    
    # draw beams
    centresx=[]
    centresy=[]
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
        centresx.append((node1_coord[0]+node2_coord[0])/2)
        centresy.append((node1_coord[1]+node2_coord[1])/2)
        fig.add_trace(go.Scatter(
            x=x,
            y=y,
            #fill='blue',
            mode="lines"
            ))
        
    #draw center points to be able to deselect beams
    fig.add_trace(go.Scatter(
        x=centresx,
        y=centresy,
        mode='markers',
        marker=dict(symbol='x')#,
        #text=np.arange(0,len(members)).astype(str)#this text only shows up when hovering
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
    quiver = create_quiver(x0,y0,fx,fy,scale=0.05,line=(dict(color='red')))
    fig.add_traces(data=quiver.data)
        
    # # draw selected point
    # if not st.session_state.selected_node == -1:
    #     sn = st.session_state.selected_node
    #     fig.add_trace(go.Scatter(x=[nodes[sn,0]],y=[nodes[sn,1]],
    #                     mode='markers',
    #                     text='you selected this one'))
    
    # Beautify plot
    #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #plt.colorbar(cpick,label="force")
    
    return fig#,ax

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
#if 'selected_member' not in st.session_state:
#    st.session_state.selected_member = []
if 'selected_nodes' not in st.session_state:
    st.session_state.selected_nodes = []

st.title("Calculating internal forces of a beam structure")

###############################################################################
# INPUTS
###############################################################################

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

fig = update_plot(internal_forces,members,nodes,f_ext)

#                 np.zeros([2*len(nodes),1])
# for i in np.arange(0,len(label_vecotr)):
#     label_vector

###############################################################################
# OUTPUTS
###############################################################################

st.write('you need to fulfill 2*n_nodes = n_members + n_supports')

with st.expander('Look at the plot', expanded=True):
    sn = plotly_events(fig)#, click_event=True)
    st.write('return value of plotly_events: ' + str(sn))
    if not sn == []:
        if sn[0]['curveNumber'] == 0:
            st.session_state.selected_nodes.append(sn[0]['pointNumber'])
            st.write('You selected node #'
                         + str(sn[0]['pointNumber']+1)
                         + '. Select another one to draw a new beam')
            st.write('current st.session_state.selected_nodes (actual node inidces, mind you): ' + str(st.session_state.selected_nodes))
            if len(st.session_state.selected_nodes) == 2:
                #st.session_state.selected_nodes[1] = sn[0]['pointNumber']
                #new_member(st.session_state.selected_nodes)
                st.session_state.selected_nodes = []
            elif len(st.session_state.selected_nodes) > 2:
                st.session_state.selected_nodes = []
            
                
                
        if sn[0]['curveNumber'] == len(members)+1:
            st.session_state.selected_member[0] = sn[0]['pointNumber']


with st.expander('Look at the Matrix. Select font size in the sidebar', expanded=True):
    st.markdown(print_equations(matrix, rhs, internal_forces,len(members),len(support),decimals,textsize))