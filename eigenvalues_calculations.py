# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 15:12:09 2022

@author: phili
"""

import numpy as np
import streamlit as st
# from printing import bmatrix

###############################################################################
# CALCULATIONS
###############################################################################

def calculate_stiffness_matrix(nodes, members, member_angles, member_lengths, A, E):
    # stiffness is calculated per member as a 4x4 matrix
    # K^(ij) = (AE/L^(ij))(...)
    def K_from_theta(theta, length, A, E):
        K_member = np.zeros((4, 4))
        c = np.cos(theta)
        s = np.sin(theta)
        c2 = c**2
        s2 = s**2
        cs = c*s
        K_member[0, 0] = c2
        K_member[1, 1] = s2
        K_member[2, 2] = c2
        K_member[3, 3] = s2
        K_member[0, 2] = -c2
        K_member[1, 3] = -s2
        K_member[2, 0] = -c2
        K_member[3, 1] = -s2
        K_member[0, 1] = cs
        K_member[1, 0] = cs
        K_member[2, 3] = cs
        K_member[3, 2] = cs
        K_member[0, 3] = -cs
        K_member[1, 2] = -cs
        K_member[2, 1] = -cs
        K_member[3, 0] = -cs
        K_member = K_member * A * E / length
        K11 = K_member[0:2,0:2]
        K12 = K_member[0:2,2:4]
        K21 = K_member[2:4,0:2]
        K22 = K_member[2:4,2:4]
        return K_member,K11,K12,K21,K22
    
    n_nodes = len(nodes)
    K = np.zeros((2*n_nodes, 2*n_nodes))
    n_members = len(members)
    members_range = np.arange(0, n_members)
    for m in members_range:
        node1 = members[m, 0]
        node2 = members[m, 1]
        theta = member_angles[m]
        length = member_lengths[m]
        K_member,K11,K12,K21,K22 = K_from_theta(theta, length, A, E)
        # st.write('member: ' + str(m) + ', ' + str(members[m]))
        # st.markdown('$' + bmatrix(K_member.round(2), ' 0 ','b') + '$')
        # st.markdown('$K11 = ' + bmatrix(K11.round(2), ' 0 ','b') + '$')
        k1 = 2*node1
        k2 = 2*node2
        K[k1:k1+2,k1:k1+2] += K11
        K[k1:k1+2,k2:k2+2] += K12
        K[k2:k2+2,k1:k1+2] += K21
        K[k2:k2+2,k2:k2+2] += K22
                
    return K

def construct_mass_matrix(node_masses):
    n_nodes = len(node_masses)
    nodes_range = np.arange(0, n_nodes)
    # n_supports = len(supports)
    # supports_range = np.arange(0, n_supports)
    # deleterows = np.empty()
    # for n in nodes_range:
    #     for s in supports_range:
    #         if n == supports[s,0]:
    #             if supports[s,1]%(pi) == 0: # horizontal movement suppressed

    M = np.zeros((2*n_nodes,2*n_nodes))
    for n in nodes_range:
        M[2*n,2*n] = node_masses[n]
        M[2*n+1,2*n+1] = node_masses[n]
        # # remove masses of nodes with supports and of those without members
        # for s in supports_range:
        #     if n == supports[s, 0]:
        #         if supports[s, 1] % (np.pi) < 1e-8:  # horizontal movement suppressed
        #             M[2*n] = 0
        #         if supports[s, 1] % (np.pi) < 1e-8:  # vertical movement suppressed
        #             M[2*n + 1] = 0
    # M = M[M > 1e-8]

    return M

# def calculate_potential_energy():
# def calculate_deformation_energy():
# def calculate_kinetic_energy(M):

def calculate_node_masses_and_member_lengths(nodes, members, rho, A, gridsize):
    n_members = len(members)
    members_range = np.arange(0, n_members)
    n_nodes = len(nodes)

    member_lengths = np.zeros(n_members)
    member_angles = np.zeros(n_members)
    for m in members_range:
        node1 = members[m, 0]
        node2 = members[m, 1]
        dx = nodes[node1][0] - nodes[node2][0]
        dy = nodes[node1][1] - nodes[node2][1]
        member_lengths[m] = np.sqrt(dx**2+dy**2)
        if dy == 0:
            div = np.inf
        else:
            div = dx/dy
        member_angles[m] = np.arctan(div)
        
    member_lengths *= gridsize

    member_masses = member_lengths * rho * A

    node_masses = np.zeros(n_nodes)
    for m in members_range:
        node1 = members[m, 0]
        node2 = members[m, 1]
        node_masses[node1] += member_masses[m]/2
        node_masses[node2] += member_masses[m]/2

    # for n in nodes_range:
    #     for m in members_range:
    #         member = members[m]
    #         if (n == member[0]) | (n == member[1]):
    #             node_masses[i] += member_masses[j]/2
    #         j += 1

    return node_masses, member_lengths, member_angles

def calculate_forces(nodes,f_ext,node_masses,gravity):
    n_nodes = len(nodes)
    nodes_range = range(n_nodes)
    force_range = range(len(f_ext))
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
    return rhs

def reduce_K(A,nodes,rods_per_node,support):
    delete_supports = True
    support_range = range(len(support))
    n_nodes = len(nodes)
    nodes_range = range(n_nodes)
    deleterows = []
    for inode in nodes_range:  # number of the node we are talking about
        k = 2*inode
        if rods_per_node[inode] < 1:
            deleterows.append(k)
            deleterows.append(k+1)
        if delete_supports:
            for s in support_range:
                if support[s][0] == inode:
                    if support[s][1] == 0:
                        deleterows.append(k)
                    if support[s][1] == 90:
                        deleterows.append(k+1)
    A = np.delete(A,deleterows,0)
    if np.shape(A)[1] > 1:
        A = np.delete(A,deleterows,1)
    return A

# def reduce_M(A):
#     keep_row = np.zeros(np.shape(A)[0],dtype=bool)
#     for i in range(np.shape(A)[0]):
#         if np.sum(A[i,:]>1e-9) > 0:
#             keep_row[i] = True
#     A = A[keep_row,:]
#     keep_column = np.zeros(np.shape(A)[1],dtype=bool)
#     for j in range(np.shape(A)[1]):
#         if np.sum(A[:,j]>1e-9) > 0:
#             keep_column[j] = True
#     A = A[:,keep_column]
#     return A