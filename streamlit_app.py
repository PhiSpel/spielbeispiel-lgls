# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:38:24 2022

@author: Philipp Spelten
"""

import numpy as np

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcol
import matplotlib.cm as cm

import plotly.express as px
#import chart_studio.plotly as py
import plotly.graph_objects as go

import streamlit as st
from streamlit_plotly_events import plotly_events

import math

###############################################################################
# CALCULATIONS
###############################################################################

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

###############################################################################
# PLOTS
###############################################################################

def update_plot(internal_forces,members,nodes):

    fig = go.Figure()
    
    fig.add_trace(go.Scatter(x=nodes[:,0],y=nodes[:,1],
                    mode='markers',
                    text=np.arange(1,len(nodes)+1)))
    
    #for n in np.arange(0,len(nodes)):
        #fig.add_trace(go.Text(nodes[n,0], nodes[n,1], str(n+1), c='red'))
    
    # Make a user-defined colormap.
    #cm1 = mcol.LinearSegmentedColormap.from_list("MyCmapName",["r","r","w","b","b"])
    # Make a normalizer that will map the time values from
    # [start_time,end_time+1] -> [0,1].
    #cnorm = mcol.Normalize(vmin=min(internal_forces),vmax=max(internal_forces))
    # Turn these into an object that can be used to map time values to colors and
    # can be passed to plt.colorbar().
    #cpick = cm.ScalarMappable(norm=cnorm,cmap=cm1)
    #cpick.set_array([])
    
    for m in np.arange(0,len(members)):
        node1 = int(members[m,0]-1)
        node2 = int(members[m,1]-1)
        A = nodes[node1,:]
        B = nodes[node2,:]
        fig.add_trace(go.Scatter(
            x=[A[0],B[0]],
            y=[A[1],B[1]],
            #fill='blue',#internal_forces[m],
            mode="lines",
            text=np.arange(1,len(members)+1)))
        #fig.add_trace(go.Text((B[0]+A[0])/2, (A[1]+B[1])/2, str(m+1)))
    
    # Beautify plot
    #ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    #plt.colorbar(cpick,label="force")
    
    #fig.show()
    
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

def print_equations(matrix, rhs, internal_forces,decimals,textsize):
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
    equation_string += r'\text{nodes} & '
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

fig = update_plot(internal_forces,members,nodes)

#                 np.zeros([2*len(nodes),1])
# for i in np.arange(0,len(label_vecotr)):
#     label_vector

###############################################################################
# OUTPUTS
###############################################################################

with st.expander('Look at the Matrix. Select font size in the sidebar', expanded=True):
    st.markdown(print_equations(matrix, rhs, internal_forces,decimals,textsize))

with st.expander('Look at the plot', expanded=False):
    #scatter = py.plot_mpl(fig, filename="my first plotly plot")
    #scatter = px.scatter(nodes[:,0],nodes[:,1])
    selected_points = plotly_events(fig)
    st.write(str(selected_points))