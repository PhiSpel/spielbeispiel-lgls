# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 13:38:24 2022

@author: Philipp Spelten
"""

from plots import update_plot
from data import check_data,update_data
from state_stuff import new_member,new_force_input,new_support_input,new_force_delete,new_support_delete,new_array,string_to_list,update_all,reset_data,setup_session_states
from printing import bmatrix,print_equations,latex_to_md
from eigenvalues_calculations import reduce_K,calculate_forces,calculate_stiffness_matrix,calculate_node_masses_and_member_lengths,construct_mass_matrix

import numpy as np

import streamlit as st
from streamlit import session_state as state
from streamlit_plotly_events import plotly_events

###############################################################################
# TO-DO
###############################################################################

# 1. Calculate Energies P, D, T
# 2. Construct Stiffness Matrix K
# 3. Two results: Deformation (P+D = min) and Eigenvalues (Mx_ddot + Kx = 0)

###############################################################################
# STREAMLIT AND STATES
###############################################################################

if True:
    st.set_page_config(layout="wide", initial_sidebar_state='collapsed')

    setup_session_states()
    
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

if True:
    debug = False  # st.sidebar.checkbox(label="show development stuff")
    
    if debug:
        decimals = st.sidebar.number_input(label="precision of print", min_value=0,max_value=5,value=2)
        # textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize',r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)
    else:
        decimals = 2
        # textsize = r'\scriptsize'
    
    if not debug:
        if 'support' not in state:
            state.support = np.array([
                [0, 0], [0, 90],[20, 90]
            ])
        if 'f_ext' not in state:
            state.f_ext = np.array([
                [5, -90,10],[12,180,10],[15,-90,15]
            ])

###############################################################################
# BUTTONS
###############################################################################

if True:
    col1, col2,col3,col4 = st.columns(4)#([1,1,1,1,0])
    with col1:
        reset = st.button('Reset data')
    if reset:
        reset_data()
    with col2:
        onlyviz = st.checkbox("Interactive mode", key='onlyviz')
    with col3:
        apply_changes = st.button('Update plot and matrix')  # ,on_click=update_members(state.removed_members,state.new_members))
    with col4:
        calculate = st.button('Update calculations')
        if calculate:
            onlyviz = False
    
    # if onlyviz:
    #     st.write('Mind the admonitions below the plot about requirements to the rod system.')

    textsize = st.sidebar.selectbox(label="font size of formula", options=[r'\normalsize', r'\small',r'\footnotesize',r'\scriptsize',r'\tiny'],index=3)
    textsize_md = latex_to_md(textsize)
    
    showzeros = st.sidebar.selectbox(label="show zeros as...", options=[' 0 ', ' '])
    
    if onlyviz:
        vectorinput = st.sidebar.checkbox(
            "Use vectors as input. (klick 'Update plot' once to initialize)")
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
                update_all(state.removed_members, state.new_members,state.removed_forces,state.removed_supports,state.new_f_ext,state.new_supports)
                state.selected_nodes = []

###############################################################################
# CALCULATIONS
###############################################################################

gravity = True
rho = 7850e-9 # kg/mm^3
gridsize = 100 # mm
d = 1 # mm
A = d*np.pi # mm^2
E = 210000 # N/mm^2

[node_masses, member_lengths, member_angles] = calculate_node_masses_and_member_lengths(state.all_nodes, state.members, rho, A, gridsize)

[rods_per_node, connected_nodes, issquare, forces_connected] = check_data(state.all_nodes, state.members,state.support,state.f_ext)

if debug:
    
    if 'matrix' not in state:
        [state.matrix, state.rhs, state.internal_forces] = update_data(state.all_nodes, state.members,state.support,state.f_ext,debug,rods_per_node,gravity,node_masses)
    if 'fig' not in state:
        [state.fig, state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
    
    if onlyviz:
        if state.onlyviz:
            [state.fig, state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
        if apply_changes:
            [state.fig, state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
    elif calculate:
        if issquare:
            [state.matrix, state.rhs, state.internal_forces] = update_data(state.all_nodes, state.members,state.support,state.f_ext,debug,rods_per_node,gravity,node_masses)
            [state.fig, state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)
        else:
            st.warning("I am having issues solving your system. Return to interactive mode to check whether your system is solvable. If you need to reset, press 'Reste data'.")
            onlyviz = True
            [state.fig, state.forcemap] = update_plot(state.internal_forces,state.members,state.all_nodes,state.f_ext,state.support,onlyviz)

###############################################################################
# OUTPUTS
###############################################################################


M = construct_mass_matrix(node_masses).round(decimals)
M = reduce_K(M,state.all_nodes,rods_per_node,state.support)

K = calculate_stiffness_matrix(state.all_nodes, state.members,member_angles,member_lengths, A, E).round(decimals)
if debug:
    st.markdown('$K_{whole} = ' + bmatrix(K, ' 0 ','b') + '$')
K = reduce_K(K,state.all_nodes,rods_per_node,state.support)

f_ext = calculate_forces(state.all_nodes,state.f_ext,node_masses,gravity)
f_ext = reduce_K(f_ext,state.all_nodes,rods_per_node,state.support)

if debug:
    st.write('node_masses: ' + str(node_masses))
    st.write('member_lengths: ' + str(member_lengths))
    st.write('member_angles: ' + str(member_angles))
    
    st.markdown('$M = ' + bmatrix(M, ' 0 ','b') + '$')
    st.markdown('$K = ' + bmatrix(K, ' 0 ','b') + '$')
    
    # st.write('dim of K: ' + str(np.shape(K)))
    # st.write('dim of M: ' + str(np.shape(M)))
    # st.write('dim of f_ext: ' + str(np.shape(f_ext)))
    
    
    # matrix, rhs, internal_forces = update_data(state.all_nodes, state.members,state.support,state.f_ext,debug,rods_per_node,gravity,node_masses,buildonly=True)
    # st.write('dim of rhs: ' + str(np.shape(rhs)))
    # st.markdown('$rhs = ' + bmatrix(rhs.round(2), ' 0 ','b') + '$')

    st.markdown('$\\text{Forces} = ' + bmatrix(f_ext.round(decimals), ' 0 ','b') + '$')

displacements = np.zeros([2*len(connected_nodes),1])
st.markdown(
    print_equations(
        K, M, f_ext, displacements, 
        len(state.members), len(state.support),
        decimals, textsize,len(connected_nodes),
        onlyviz, showzeros,connected_nodes,rods_per_node))

left_side = K#[2:,2:]
right_side = -np.dot(M,f_ext)#[2:]
displacements = np.linalg.solve(left_side,right_side)
st.markdown('$\\text{displacements}: ' + bmatrix(displacements.round(4), ' 0 ','b') + '\\text{mm}$')

#m_to_half = np.linalg.sqrtm(M)
#m_to_minus_half = np.linalg.inv(m_to_half)
#m_to_minus_half = np.linalg.matrix_power(M,-0.5)
m_to_minus_half = np.power(np.sqrt(M),-1)
m_to_minus_half[m_to_minus_half == np.inf] = 0
ev_problem = np.dot(np.dot(m_to_minus_half,K),m_to_minus_half)

w, v = np.linalg.eig(ev_problem)

#st.write(str(w))
#st.write(str(v))
st.markdown('$w = ' + bmatrix([w.round(4)], ' 0 ','b') + '$')
st.markdown('$v = ' + bmatrix(v.round(4), ' 0 ','b') + '$')


###############################################################################
# OLD OUTPUTS
###############################################################################

# loading fig from cache
# fig = state.fig

# if onlyviz:
#     sn = plotly_events(fig)
#     if debug:
#         st.write('return value of plotly_events: ' + str(sn))
#     if not sn == []:

#         if sn[0]['curveNumber'] == len(state.members)+0:
#             state.selected_nodes.append(sn[0]['pointNumber'])
#             if len(state.selected_nodes) == 1:
#                 st.markdown(textsize_md + 'You selected node #'
#                             + str(sn[0]['pointNumber'])
#                             + '. Select another one to draw a new rod. '
#                             + r''' </font>''', unsafe_allow_html=True)
#             elif len(state.selected_nodes) == 2:
#                 if not state.selected_nodes[1] == state.selected_nodes[0]:
#                     new_member(state.selected_nodes)
#                     st.markdown(textsize_md + 'You selected nodes #'
#                                 + str(state.new_members[-1])
#                          + ". See below which changes will be applied by 'update plot'. "
#                          + r''' </font>''', unsafe_allow_html=True)
#                     state.selected_nodes = []
#                 else:
#                     st.markdown(textsize_md + 'You selected node #'
#                                 + str(sn[0]['pointNumber'])
#                                 + '. Twice. Select a different one. '
#                                 + r''' </font>''', unsafe_allow_html=True)
#             elif len(state.selected_nodes) > 2:
#                 state.selected_nodes = []
#                 st.markdown(textsize_md + 'You selected more than 2 nodes. Must be a bug. Cleared selected nodes.'
#                             + r''' </font>''', unsafe_allow_html=True)
#             if debug:
#                 nodes_string = 'Current storage in state.selected_nodes: ' + \
#                     str(state.selected_nodes)
#                 if len(state.selected_nodes) == 2:
#                     if state.selected_nodes[1] == state.selected_nodes[0]:
#                         nodes_string += '. Will not write to state.new_members, b/c you selected the same node twice.'
#                     else:
#                         new_member(state.selected_nodes)
#                         nodes_string += '. Wrote nodes to state.new_members.'
#                     state.selected_nodes = []
#                 elif len(state.selected_nodes) > 2:
#                     state.selected_nodes = []
#                     nodes_string += '. Must be a bug. Cleared selected nodes.'
#                 st.write(nodes_string)

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

        # if sn[0]['curveNumber'] == len(state.members)+3:
        #     if sn[0]['pointNumber'] not in state.removed_members:
        #         state.removed_members.append(sn[0]['pointNumber'])

    # checking if we will get a square matrix
#     if forces_connected == False:
#         st.warning(
#             'Not all forces are connected. You may solve the system anyway.')

#     st.markdown(textsize_md
#                 + r"'Update plot' will remove **members** #" + \
#                     str(state.removed_members)
#                 + ' and add members on ' + str(state.new_members) + '. '
#                 + r'You remove **supports** #' + str(state.removed_supports)
#                 + 'and add supports on ' + str(np.round(state.new_supports, 2)) + '. '
#                 + r'You remove **forces** #' + str(state.removed_forces)
#                 + ' and add forces on ' + str(np.round(state.new_f_ext, 2)) + '. '
#                 + r''' </font>''', unsafe_allow_html=True)

#     status_string = (textsize_md + r'''You need to fulfill 
#         $$ 2 \cdot n_\text{nodes} = n_\text{members} + n_\text{supports} $$
#         to get a square matrix. ''' + 'You have '
#                      + str(len(connected_nodes)) + ' nodes, '
#         + str(len(state.members)) + ' members and '
#         + str(len(state.support)) + ' supports. ')
#     if issquare:
#         status_string += 'Currently, you get a square matrix.'
#     else:
#         if (2*len(connected_nodes) > (len(state.members) + len(state.support))):
#             status_string += 'You need to add more members or supports.'
#         else:
#             status_string += 'You need to remove members or supports or add more nodes.'
#     st.markdown(status_string + r''' </font>''', unsafe_allow_html=True)

# else:
#     st.plotly_chart(fig)

# if onlyviz:
#     buildonly = True
#     matrix, rhs, internal_forces = update_data(state.all_nodes, state.members,state.support,state.f_ext,debug,rods_per_node,buildonly)
#     st.markdown(
#         print_equations(
#             matrix, rhs, internal_forces,
#             len(state.members), len(state.support),
#             decimals, textsize,len(connected_nodes),
#             onlyviz, showzeros,connected_nodes,rods_per_node))
# else:
#     st.markdown(
#         print_equations(
#             state.matrix, state.rhs, state.internal_forces,
#             len(state.members), len(state.support),
#             decimals, textsize,len(connected_nodes),
#             onlyviz, showzeros,connected_nodes,rods_per_node))

# truss1 = Image.open('images/truss1.png')
content = """<p><a href='#' id='Link 1'>First link</a></p>
<p><a href='#' id='Link 2'>Second link</a></p>
<a href='#' id='Image 1'><img width='20%' src='https://images.unsplash.com/photo-1565130838609-c3a86655db61?w=200'></a>
<a href='#' id='truss1'><img width='20%' src='https://github.com/PhiSpel/spielbeispiel-lgls/blob/master/images/truss1.png?raw=true'></a>
"""
