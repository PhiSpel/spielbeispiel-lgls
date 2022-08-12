# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:51:06 2022

@author: phili
"""

###############################################################################
# DATA WORK
###############################################################################

import numpy as np
from streamlit import session_state as state
import streamlit as st
from printing import print_equations


def check_data(nodes, members,support,f_ext):
    n_members = len(members)
    members_range = np.arange(0, n_members)
    n_nodes = len(nodes)
    nodes_range = np.arange(0, n_nodes)
    nsupport = len(support)
    support_range = np.arange(0, nsupport)
    nf_ext = len(f_ext)
    force_range = np.arange(0, nf_ext)
    rods_per_node = np.zeros(n_nodes)
    for m in members_range:
        node1 = members[m, 0]
        node2 = members[m, 1]
        rods_per_node[node1] += 1
        rods_per_node[node2] += 1
    for s in support_range:
        inode = int(support[s, 0])
        rods_per_node[inode] += 1
    for f in force_range:
        inode = int(f_ext[f, 0])
        rods_per_node[inode] += 1
    deletenodes = []
    for inode in nodes_range:
        if rods_per_node[inode] < 2:
            deletenodes.append(inode)
    connected_nodes = np.delete(nodes,deletenodes,0)
    # st.write(str(connected_nodes))
    forces_connected = True
    for f in force_range:
        if int(f_ext[f, 0]) in deletenodes:
            forces_connected = False

    issquare = 2*len(connected_nodes) == (len(state.members) + \
                     len(state.support))

    return rods_per_node, connected_nodes, issquare, forces_connected

# maybe use @st.cache()?
def update_data(nodes, members,support,f_ext,debug,rods_per_node,gravity,node_masses,buildonly=False):

    n_members = len(members)
    members_range = np.arange(0, n_members)
    n_nodes = len(nodes)
    nodes_range = np.arange(0, n_nodes)
    nsupport = len(support)
    support_range = np.arange(0, nsupport)
    nf_ext = len(f_ext)
    force_range = np.arange(0, nf_ext)

    members = np.append(members, np.zeros([n_members, 3]),axis = 1)
    # compute angles
    for m in members_range:
        node1 = int(members[m, 0])
        node2 = int(members[m, 1])
        stab = nodes[node2, :]-nodes[node1,:]
        direction = np.sign(stab[1])
        if direction == 0:
            direction = 1
        alpha = np.arccos(stab[0]/np.linalg.norm(stab))*direction
        members[m, 2] = np.cos(alpha) # along x
        members[m, 3] = np.sin(alpha) # along y
        members[m, 4] = alpha

    # compute matrix coefficients
    ks1 = np.zeros([n_nodes, n_members])
    ks2 = np.zeros([n_nodes, n_members])
    matrix = np.zeros([2*n_nodes, n_members])
    for k in nodes_range:
        # compute x and y forces per member
        for m in members_range:
            node1 = int(members[m, 0])
            node2 = int(members[m, 1])
            ks1[k,m] = int(k == node1)
            ks2[k,m] = -int(k == node2)
        ks = ks1+ks2
        # forces along x
        matrix[2*k, :] = np.multiply(ks[k,:],members[:,2])
        # forces along y
        matrix[2*k+1, :] = np.multiply(ks[k,:],members[:,3])

    # add support forces to matrix
    if debug:
        st.write('supports: ' + str(support))
    matrix = np.append(matrix, np.zeros([len(matrix), len(support)]),axis = 1)
    for s in support_range:
        inode = int(support[s, 0])
        angle = np.radians(support[s, 1])
        if debug:
            st.write('writing support #' + str(s) + \
                     ' at node #' + str(inode) + ' to matrix')
        matrix[2*inode, n_members+s] = np.cos(angle) # forces along x
        matrix[2*inode+1, n_members+s] = np.sin(angle) # forces along y

    # compute right hand side
    rhs = np.zeros([2*n_nodes, 1])
    for force in force_range:
        inode = int(f_ext[force, 0])
        angle = np.radians(f_ext[force, 1])
        newtons = f_ext[force, 2]
        # external force along x
        rhs[2*inode, 0] = newtons*np.cos(angle)
        # external force along y
        rhs[2*inode+1, 0] = newtons*np.sin(angle)
        
    if gravity:
        for n in nodes_range:
            rhs[2*n + 1] += -9.81*node_masses[n]
            # st.write('node: ' + str(n) + ', weight: ' + str(node_masses[n]))

    # delete obsolete rows
    deleterows = []
    for inode in nodes_range:  # number of the node we are talking about
        k = 2*inode
        if rods_per_node[inode] <= 1:
            deleterows.append(k)
            deleterows.append(k+1)
    if debug:
        # connected_nodes=np.delete(nodes,deletenodes,0)
        st.write('nodes: ' + str(nodes) + ' n_nodes: ' + str(n_nodes))
        st.markdown(print_equations(matrix, rhs, [], len(members),len(state.support),1,'\scriptsize'))
        #st.write('connected nodes: ' + str(connected_nodes))
        st.write('delete rows: ' + str(deleterows))
    matrix = np.delete(matrix, deleterows,0)
    rhs = np.delete(rhs, deleterows,0)    
    # st.write(str(matrix))

    if debug:
        st.markdown(print_equations(matrix, rhs, [], len(members),len(state.support),1,'\scriptsize'))

    # solve for unknowns
    if buildonly:
        internal_forces = []
    else:
        internal_forces = np.linalg.solve(matrix, rhs)

    return matrix, rhs, internal_forces