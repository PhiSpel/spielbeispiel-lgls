# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:57:42 2022

@author: phili
"""

import numpy as np

from streamlit import session_state as state

###############################################################################
# PRINTING OUTPUT
###############################################################################

def bmatrix(a, showzeros,matrixtype=''):
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

def print_equations(K, M, f_ext, displacements, n_rods,n_bcs,decimals,textsize,n_nodes,onlyviz,showzeros,connected_nodes,rods_per_node):
    from eigenvalues_calculations import reduce_K
    
    label_vector = []
    k = 0
    for node in state.all_nodes:
        label_vector.append([r' \text{'
                             # + 'N.' + str(int(np.floor(k/2))) + ' '
                             + str(node)])
        label_vector[k][0] += ', x'
        label_vector[k][0] += '}'
        k += 1
        label_vector.append([r' \text{'
                             # + 'N.' + str(int(np.floor(k/2))) + r' '
                             + str(node)])
        label_vector[k][0] += ', y'
        label_vector[k][0] += '}'
        k += 1
    label_vector = np.array(label_vector)
    label_vector = reduce_K(label_vector,state.all_nodes,rods_per_node,state.support)
    
    def rounding(a):
        a = np.round(a, decimals)
        a[a == 0] = 0
        return a
    
    M_f             = -np.dot(M,f_ext)
    M_f             = rounding(M_f)
    M               = rounding(M)
    f_ext           = rounding(f_ext)
    K               = rounding(K)
    displacements   = rounding(displacements)
    
    equation_string = r'$'
    equation_string += textsize
    equation_string += r' \begin{matrix} '
    equation_string += r'\text{nodes} & \text{stiffness matrix}'
    # equation_string += bmatrix(label_vector,'h')
    equation_string += r' & \text{displacements} & -M\cdot f_{ext} & & \\'
    equation_string += bmatrix(label_vector, showzeros,'v')
    equation_string += ' & '
    equation_string += bmatrix(K, showzeros,'b')
    equation_string += ' & '
    equation_string += '\cdot '
    equation_string += bmatrix(displacements, showzeros,'b')
    equation_string += ' & '
    equation_string += '='
    equation_string += bmatrix(-M_f, showzeros,'b')
    equation_string += ' & '
    equation_string += '='
    equation_string += bmatrix(-M, showzeros,'b')
    equation_string += ' & '
    equation_string += '\cdot '
    equation_string += bmatrix(f_ext, showzeros,'b')
    equation_string += '\end{matrix}$'
    return equation_string

# def print_equations(matrix, rhs, internal_forces, n_rods,n_bcs,decimals,textsize,n_nodes,onlyviz,showzeros,connected_nodes):
#     label_vector = []
#     k = 0
#     for node in connected_nodes:
#         label_vector.append([r' \text{'
#                              # + 'N.' + str(int(np.floor(k/2))) + ' '
#                              + str(node)])
#         label_vector[k][0] += ', x'
#         label_vector[k][0] += '}'
#         k += 1
#         label_vector.append([r' \text{'
#                              # + 'N.' + str(int(np.floor(k/2))) + r' '
#                              + str(node)])
#         label_vector[k][0] += ', y'
#         label_vector[k][0] += '}'
#         k += 1
#     label_vector = np.array(label_vector)
#     matrix = np.round(matrix, decimals)
#     matrix[matrix == 0] = 0
#     rhs = np.round(rhs, decimals)
#     rhs[rhs == 0] = 0
#     internal_forces = np.round(internal_forces, decimals)
#     internal_forces[internal_forces == 0] = 0
#     equation_string = r'$'
#     equation_string += textsize
#     equation_string += r' \begin{matrix} '
#     equation_string += r'\text{' + str(n_nodes) + r' nodes} & \text{stiffness matrix}'
#     # equation_string += bmatrix(label_vector,'h')
#     equation_string += r' & \text{displacements} & \text{forces}\\'
#     equation_string += bmatrix(label_vector, showzeros,'v')
#     equation_string += ' & '
#     equation_string += bmatrix(matrix, showzeros,'b')
#     equation_string += ' & '
#     equation_string += '\cdot '
#     if (onlyviz) | (internal_forces == []):
#         equation_string += '?'
#     else:
#         equation_string += bmatrix(internal_forces, showzeros,'b')
#     equation_string += ' & '
#     equation_string += '='
#     equation_string += bmatrix(rhs, showzeros,'b')
#     equation_string += '\end{matrix}$'
#     return equation_string

def latex_to_md(textsize):
    string = r'''<font size="'''
    if textsize == r'\tiny':
        string += '1'
    elif textsize == r'\scriptsize':
        string += '2'
    elif (textsize == r'\footnotesize') | (textsize == r'\small'):
        string += '3'
    elif textsize == r'\normalsize':
        string += '4'
    string += r'''"> '''
    return string