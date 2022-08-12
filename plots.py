# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 14:48:07 2022

@author: phili
"""

###############################################################################
# PLOTS
###############################################################################

import numpy as np

import plotly.express as px
# import chart_studio.plotly as py
import plotly.graph_objects as go
import plotly
# import plotly.colors.sequential as seq_color
from plotly.figure_factory import create_quiver
# Identical to Adam's answer
import plotly.colors

def get_continuous_color(colorscale, intermed):
    '''# Plotly continuous colorscales assign colors to the range [0, 1]. This function computes the intermediate
    # color for any value in that range.

    # Plotly doesn't make the colorscales directly accessible in a common format.
    # Some are ready to use:

    #     colorscale = plotly.colors.PLOTLY_SCALES["Greens"]

    # Others are just swatches that need to be constructed into a colorscale:

    #     viridis_colors, scale = plotly.colors.convert_colors_to_same_type(plotly.colors.sequential.Viridis)
    #     colorscale = plotly.colors.make_colorscale(viridis_colors, scale=scale)

    # :param colorscale: A plotly continuous colorscale defined with RGB string colors.
    # :param intermed: value in the range [0, 1]
    # :return: color in rgb string format
    # :rtype: str
    '''
    if len(colorscale) < 1:
        raise ValueError("colorscale must have at least one color")

    def hex_to_rgb(c): return "rgb" + str(ImageColor.getcolor(c, "RGB"))

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

def update_plot(internal_forces, members,nodes,f_ext,support,onlyviz):

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

    max_force = max(abs(internal_forces))[0]
    normed_force = (internal_forces/max_force + 1)/2

    for m in np.arange(0, len(members)):
        node1_index = int(members[m, 0])
        node2_index = int(members[m, 1])
        node1_coord = nodes[node1_index, :]
        node2_coord = nodes[node2_index, :]
        x = []
        y = []
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
                    color=get_color('rdbu', normed_force[m][0])
                ),
                marker=dict(
                    cmax=max_force,
                    cmin=-max_force,
                    colorbar=dict(
                        title="Force intensity",
                        orientation='h'  # this does nothing! :/
                    ),
                    colorscale="rdbu"
                )
            ))

    # create a heatmap to be added below
    forcemap = px.imshow([np.linspace(-max_force, max_force,100)],
                         color_continuous_scale='rdbu')
    forcemap.layout.coloraxis.showscale = False

    # draw nodes last, so they can be selected
    fig.add_trace(go.Scatter(x=nodes[:, 0],y=nodes[:,1],
                             mode='markers',
                    text=np.arange(0, len(nodes)),
                             name='nodes'
                             ))

    # calculate coordinates of forces
    x0 = []
    y0 = []
    fx = []
    fy = []
    for f in np.arange(0, len(f_ext)):
        node_index = int(f_ext[f, 0])
        x0.append(nodes[node_index, 0])
        y0.append(nodes[node_index, 1])
        angle = np.radians(f_ext[f, 1])
        newtons = f_ext[f, 2]
        # external force along x
        fx.append(newtons*np.cos(angle))
        # external force along y
        fy.append(newtons*np.sin(angle))
    # draw forces
    quiver = create_quiver(x0, y0,fx,fy,scale=0.05,line=(dict(color='red')),name='Forces')
    fig.add_traces(data=quiver.data)

    # draw supports
    support_length = 5
    x0 = []
    y0 = []
    sx = []
    sy = []
    for s in np.arange(0, len(support)):
        node_index = int(support[s, 0])
        x0.append(nodes[node_index, 0])
        y0.append(nodes[node_index, 1])
        angle = np.radians(support[s, 1])
        # external force along x
        sx.append(support_length*np.cos(angle))
        # external force along y
        sy.append(support_length*np.sin(angle))
    # draw forces
    quiver = create_quiver(x0, y0,sx,sy,
                           scale=0.05,
                           line=(dict(color='green')),
                           name='Supports')
    fig.add_traces(data=quiver.data)

    # draw center points to be able to deselect rods
    centresx = []
    centresy = []
    member_names = []
    for m in np.arange(0, len(members)):
        node1_index = int(members[m, 0])
        node2_index = int(members[m, 1])
        node1_coord = nodes[node1_index, :]
        node2_coord = nodes[node2_index, :]
        centresx.append((node1_coord[0]+node2_coord[0])/2)
        centresy.append((node1_coord[1]+node2_coord[1])/2)
        if onlyviz:
            member_names.append('deselect member #' + str(m))
        else:
            member_names.append(
                'member #' + str(m) + '; current force: ' + str(round(internal_forces[m][0])))
    if onlyviz:
        fig.add_trace(go.Scatter(
            x=centresx,
            y=centresy,
            mode='markers',
            marker=dict(symbol='x', color='orange'),
            name='members for deselection',
            text=member_names
        ))
    else:
        fig.add_trace(go.Scatter(
            x=centresx,
            y=centresy,
            mode='markers',
            marker=dict(symbol='circle', color='orange'),
            name='',
            text=member_names
        ))

    return fig, forcemap