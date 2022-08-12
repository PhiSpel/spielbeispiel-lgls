# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:53:22 2022

@author: phili
"""

###############################################################################
# NEW STUFF
################################################################################

import streamlit as st
from streamlit import session_state as state
import numpy as np

def new_member(new_nodes):
    node1 = new_nodes[0]
    node2 = new_nodes[1]
    if state.new_members == [[]]:
        state.new_members = [[node1, node2]]
    else:
        state.new_members = np.append(state.new_members, [[node1,node2]],axis=0)
    return

def new_force_input():
    if not state.new_f_ext_str == '':
        new_f_ext_list = string_to_list(state.new_f_ext_str)
        # new_f_ext_list[1] = math.radians(new_f_ext_list[1])
        if state.new_f_ext == [[]]:
            state.new_f_ext = [new_f_ext_list]
        else:
            np.append(state.new_f_ext, [new_f_ext_list],axis=0)
        state.new_f_ext_str = ''
    return

def new_support_input():
    if not state.new_support_str == '':
        new_support_list = string_to_list(state.new_support_str)
        # new_support_list[1] = math.radians(new_support_list[1])
        if state.new_supports == [[]]:
            state.new_supports = [new_support_list]
        else:
            if new_support_list not in state.new_supports:
                state.new_supports = np.append(state.new_supports, [new_support_list],axis=0)
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

def update_all(removed_members, new_members,removed_forces,removed_supports,new_f_ext,new_supports):
    if not removed_members == []:
        state.members = np.delete(state.members, removed_members,axis=0)
        state.removed_members = []
    if not new_members == [[]]:
        state.members = np.append(state.members, new_members,axis=0)
        state.new_members = [[]]

    if not removed_forces == []:
        state.f_ext = np.delete(state.f_ext, removed_forces,axis=0)
        state.removed_forces = []
    if not new_f_ext == [[]]:
        # if debug:
        #     st.write(str(new_f_ext))
        state.f_ext = np.append(state.f_ext, new_f_ext,axis=0)
        state.new_f_ext = [[]]

    if not removed_supports == []:
        state.support = np.delete(state.support, removed_supports,axis=0)
        state.removed_supports = []
    if not new_supports == [[]]:
        state.support = np.append(state.support, new_supports,axis=0)
        state.new_supports = [[]]

    return

def reset_data():
    state.members = np.array([
        [0, 6],
        [0, 5],
        [6, 5],
        [5, 10],
        [6, 10],
        [6, 12],
        [12, 10],
        [10, 15],
        [10, 16],
        [16, 15],
        [12, 16],
        [15, 20],
        [16, 20]
    ])
    state.support = np.array([
        [0, 0], [0, 90],[20, 90]
    ])
    state.f_ext = np.array([
        [5, -90,10],[12,180,10],[15,-90,15]
    ])
    return

def setup_session_states():

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
    
    if 'all_nodes' not in state:
        l = 5  # st.sidebar.number_input(label='length of plot',min_value=2,max_value=10,value=5)
        h = 4  # st.sidebar.number_input(label='height of plot',min_value=1,max_value=10,value=4)
        all_nodes = []
        for x in range(l+1):
            for y in range(h+1):
                all_nodes.append([x, y])
        state.all_nodes = np.array(all_nodes)
    
    if 'members' not in state:
        state.members = np.array([
            [0, 6],
            [0, 5],
            [6, 5],
            [5, 10],
            [6, 10],
            [6, 12],
            [12, 10],
            [10, 15],
            [10, 16],
            [16, 15],
            [12, 16],
            [15, 20],
            [16, 20]
        ])
    return