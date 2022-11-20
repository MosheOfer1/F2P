import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
import numpy as np
from matplotlib.widgets import Slider

MARKER_SIZE = 16
MARKER_SIZE_SMALL = 1
MARKERS_GAP = 0.2
LINE_WIDTH = 3
LINE_WIDTH_SMALL = 1
FONT_SIZE = 12
FONT_SIZE_SMALL = 5
LEGEND_FONT_SIZE = 8
LEGEND_FONT_SIZE_SMALL = 5
UNIFORM_CHAIN_MIG_COST = 600
DELTA_JUMPS = 0.01


# Set the parameters of the plot (sizes of fonts, legend, ticks etc.).
# mfc='none' makes the markers empty.
def set_plt_params(size='large'):
    if size == 'large':
        matplotlib.rcParams.update_graph2({'font.size': FONT_SIZE,
                                           'legend.fontsize': LEGEND_FONT_SIZE,
                                           'xtick.labelsize': FONT_SIZE,
                                           'ytick.labelsize': FONT_SIZE,
                                           'axes.labelsize': FONT_SIZE,
                                           'axes.titlesize': FONT_SIZE,
                                           'lines.marker': 'x',
                                           'lines.markeredgecolor': 'black'})
    else:
        matplotlib.rcParams.update_graph2({'font.size': FONT_SIZE_SMALL,
                                           'legend.fontsize': LEGEND_FONT_SIZE_SMALL,
                                           'xtick.labelsize': FONT_SIZE_SMALL,
                                           'ytick.labelsize': FONT_SIZE_SMALL,
                                           'axes.labelsize': FONT_SIZE_SMALL,
                                           'axes.titlesize': FONT_SIZE_SMALL, })
#
#
# def pre_calculate_stages(exp, cnt_size):
#     stage = [1]
#     for i in range(1, 2 ** exp):
#         stage.append(stage[i - 1] + 2 ** (cnt_size - exp))
#     return stage


def absolute_resolution(x_arr):
    y_arr = [1]
    for c in range(1, len(x_arr)):
        y_arr.append(x_arr[c] - x_arr[c - 1])
    return y_arr


def relative_resolution(x_arr, y_abs):
    y_arr = [1]
    for c in range(0, len(x_arr) - 1):
        y_arr.append(y_abs[c] / x_arr[c + 1])
    return y_arr


def dynamic_sead(cnt_size):
    # stage = pre_calculate_stages(cnt_size, cnt_size)
    xs_dynamic = []
    for i in range((2 ** cnt_size) - 1):
        exp_size = 0
        bin_cntr = np.binary_repr(i, cnt_size)
        for j in bin_cntr:
            if j == '0':
                break
            else:
                exp_size += 1
        mantissa = int(bin_cntr[0:cnt_size - exp_size], 2)
        xs_dynamic.append((mantissa * (2 ** exp_size))) #+ stage[exp_size])
    xs_dynamic = list(set(xs_dynamic))
    ys_dynamic_abs = absolute_resolution(xs_dynamic)
    ys_dynamic_relative = relative_resolution(xs_dynamic, ys_dynamic_abs)
    return xs_dynamic, ys_dynamic_abs, ys_dynamic_relative


def static_sead(cnt_size, exp_size):
    xs_static = []
    for e in range(2 ** exp_size):
        for m in range(2 ** (cnt_size - exp_size)):
            xs_static.append(m * 2 ** e)
    xs_static = list(set(xs_static))
    xs_static.sort()
    ys_static_abs = absolute_resolution(xs_static)
    ys_static_relative = relative_resolution(xs_static, ys_static_abs)
    return xs_static, ys_static_abs, ys_static_relative


# function that computes the shared estimators and D by the CEDAR's formula
def cedar(delta, max_val):
    shared_estimators = []
    different = []
    shared_estimators.append(0)
    i = 0
    while shared_estimators[i] < max_val:
        # the formula
        different.append((1 + 2 * delta ** 2 * shared_estimators[i]) / (1 - delta ** 2))
        i += 1
        shared_estimators.append(shared_estimators[i - 1] + different[i - 1])
    shared_estimators.pop()
    return shared_estimators, different


def ideal_exp_size(counted_num, cnt_size):
    for ideal_exp in range(1, cnt_size):
        if ((2 ** (cnt_size - ideal_exp)) - 1) * 2 ** ((2 ** ideal_exp) - 1) >= counted_num:
            return ideal_exp
    return cnt_size


# brute force method for finding the min delta in jumps of 0.01 each time by recursion
def min_delta(delta, max_val, cnt_size):
    shared_estimators = []
    different = []
    shared_estimators.append(0)
    i = 0
    while shared_estimators[i] < max_val:
        different.append((1 + 2 * delta ** 2 * shared_estimators[i]) / (1 - delta ** 2))
        i += 1
        shared_estimators.append(shared_estimators[i - 1] + different[i - 1])
    if i >= 2 ** cnt_size:
        return delta + DELTA_JUMPS
    return min_delta(delta - DELTA_JUMPS, max_val, cnt_size)


def update_byte(val):
    update_graph1(8, val)


# plot the static and dynamic sead to one graph
def plot_graphs(cnt_size, exp_size):
    figure, axis = plt.subplots(2, 2)
    Sliders.axis = axis
    byte_slider = Slider(plt.axes([0.53, 0.95, 0.1, 0.04]), 'Exp', 3, 6, 3)
    byte_slider.valstep = 1
    byte_slider.on_changed(update_byte)
    update_graph1(cnt_size, exp_size)
    plt.tight_layout()
    plt.show()


def update_graph1(cnt_size, exp_size):
    axis = Sliders.axis
    xs_static, ys_static_abs, ys_static_relative = static_sead(cnt_size, exp_size)
    # xs_dynamic, ys_dynamic_abs, ys_dynamic_relative = dynamic_sead(cnt_size, exp_size)
    xs_dynamic, ys_dynamic_abs, ys_dynamic_relative = dynamic_sead(cnt_size)

    axis[0, 0].clear()
    # axis[0, 0] = plt.gca()
    # axis[0, 0].grid(True)
    axis[0, 0].plot(xs_static, ys_static_abs, marker='^', markevery=MARKERS_GAP, linestyle='dashed')
    axis[0, 0].set_title('Static SEAD')
    axis[0, 0].set_xlabel('Real Value')
    axis[0, 0].set_ylabel('Absolute Resolution')
    axis[0, 0].set_xlim([0, xs_static[-1]])
    axis[0, 0].set_ylim(0)

    axis[1, 0].clear()
    axis[1, 0].plot(xs_static, ys_static_relative, 'C3', marker='v', markevery=MARKERS_GAP, linestyle='dashed',
                    mfc='none')
    axis[1, 0].set_xlabel('Real Value')
    axis[1, 0].set_ylabel('Relative Resolution')
    axis[1, 0].set_ylim([0, 0.15])
    axis[1, 0].set_xlim([0, xs_static[-1]])

    axis[0, 1].clear()
    axis[0, 1].plot(xs_dynamic, ys_dynamic_abs, marker='^', markevery=MARKERS_GAP, linestyle='dashed')
    axis[0, 1].set_title('Dynamic SEAD')
    axis[0, 1].set_xlabel('Real Value')
    axis[0, 1].set_ylabel('Absolute Resolution')
    axis[0, 1].set_xlim([0, xs_dynamic[-1]])
    axis[0, 1].set_ylim(0)

    axis[1, 1].clear()
    axis[1, 1].plot(xs_dynamic, ys_dynamic_relative, 'C3', marker='v', markevery=MARKERS_GAP, linestyle='dashed',
                    mfc='none')
    axis[1, 1].set_xlabel('Real Value')
    axis[1, 1].set_ylabel('Relative Resolution')
    axis[1, 1].set_ylim([0, 0.15])
    axis[1, 1].set_xlim([0, xs_dynamic[-1]])


class Sliders:
    axis = None
    axs = None
    CURRENT_DELTA = 0.1
    CURRENT_MAX_VAL = 2048


def delta_changed(val):
    Sliders.CURRENT_DELTA = val
    update_graph2(val, Sliders.CURRENT_MAX_VAL)


def max_value_changed(val):
    update_graph2(Sliders.CURRENT_DELTA, val)


def update_graph2(delta, max_val):
    axis = Sliders.axs
    xs_cedar_static, ys_cedar_static = cedar(delta, max_val)
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    # find fair conditions to compare between the different methods
    cnt_size_needed = math.ceil(math.log2(len(xs_cedar_static)))
    exp_size_needed = ideal_exp_size(xs_cedar_static[-1], cnt_size_needed)

    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size_needed, exp_size_needed)
    xs_sead_dynamic, ys_sead_dynamic_abs, ys_sead_dynamic_relative = dynamic_sead(cnt_size_needed)

    axis[0].clear()
    axis[0].plot(xs_cedar_static, ys_cedar_static, "-b", label="Static cedar", marker='^', markevery=MARKERS_GAP,
                 linestyle='dashdot')
    axis[0].plot(xs_sead_static, ys_sead_static_abs, "-r", label="Static SEAD", marker='^', markevery=MARKERS_GAP,
                 linestyle='dashed')
    axis[0].plot(xs_sead_dynamic, ys_sead_dynamic_abs, "-g", label="Dynamic SEAD", marker='v', markevery=MARKERS_GAP,
                 mfc='none')
    axis[0].legend(loc="upper left")
    axis[0].set_title(str(exp_size_needed) + " Exp bites and " + str(cnt_size_needed) +
                      " bits need in order to count up to " + str(math.ceil(xs_cedar_static[-1])) + " ")
    axis[0].set_xlabel('Real Value')
    axis[0].set_ylabel('Absolute Resolution')
    axis[0].set_xlim([0, xs_cedar_static[-1]])
    axis[0].set_ylim(0, 120)

    axis[1].clear()
    axis[1].plot(xs_cedar_static, ys_cedar_relative, "-b", label="Static cedar", marker='^', markevery=MARKERS_GAP,
                 linestyle='dashdot')
    axis[1].plot(xs_sead_static, ys_sead_static_relative, "-r", label="Static SEAD", marker='^', markevery=MARKERS_GAP,
                 linestyle='dashed')
    axis[1].plot(xs_sead_dynamic, ys_sead_dynamic_relative, "-g", label="Dynamic SEAD", marker='v',
                 markevery=MARKERS_GAP, mfc='none')
    axis[1].set_xlim([10, xs_cedar_static[-1]])
    axis[1].set_ylim(0, 0.18)
    axis[1].set_xlabel('Real Value')
    axis[1].set_ylabel('Relative Resolution')


# plot the CEDAR and the SEAD at the same graph
def plot_graph2(delta, max_val):
    fig, axs = plt.subplots(2)
    Sliders.axs = axs
    delta_slider = Slider(plt.axes([0.1, 0.95, 0.1, 0.04]), 'Delta', DELTA_JUMPS, 0.2, 0.1)
    delta_slider.valstep = DELTA_JUMPS
    delta_slider.on_changed(delta_changed)
    max_value = Slider(plt.axes([0.32, 0.95, 0.1, 0.04]), 'Max Value', 255, 5000, 2048)
    max_value.valstep = 20
    max_value.on_changed(max_value_changed)
    update_graph2(delta, max_val)
    plt.tight_layout()
    plt.show()


def f2p_VS_rest(axis, cnt_size, x_counter, y_counter, f2p_axis_abs, f2p_axis_rel):
    xs_f2p, ys_f2p_abs, ys_f2p_relative = axis
    # for fare comparison with CEDAR we need to check the min delta
    min_del = min_delta(0.1, xs_f2p[-1], cnt_size)
    # for fare comparison with CEDAR we need to check the min exp
    min_exp = ideal_exp_size(xs_f2p[-1], cnt_size)
    # compare to SEAD method and CEDAR to F2P
    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size, min_exp)
    xs_cedar_static, ys_cedar_static = cedar(min_del, xs_f2p[-1])
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    f2p_axis_abs[x_counter, y_counter].set_title(
        str(cnt_size) + " bits, " + str(round(min_del, 2)) + " Delta, " + str(min_exp) + " exp, " + "max " + str(
            xs_f2p[-1]))
    f2p_axis_abs[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_static, "-b", label="Static cedar")
    f2p_axis_abs[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_abs, "-r", label="Static SEAD")
    f2p_axis_abs[x_counter, y_counter].plot(xs_f2p, ys_f2p_abs, "-gD", label="F2P", markevery=0.2)
    f2p_axis_abs[x_counter, y_counter].set_xlim([0, xs_cedar_static[-1]])
    f2p_axis_abs[x_counter, y_counter].set_ylim([0, ys_f2p_abs[-1] * 5])
    plt.tight_layout()

    f2p_axis_rel[x_counter, y_counter].set_title(
        str(cnt_size) + " bits, " + str(round(min_del, 2)) + " Delta, " + str(min_exp) + " exp, " + "max " + str(
            xs_f2p[-1]))
    f2p_axis_rel[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_relative, "-b", label="Static cedar")
    f2p_axis_rel[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_relative, "-r", label="Static SEAD")
    f2p_axis_rel[x_counter, y_counter].plot(xs_f2p, ys_f2p_relative, "-gD", label="F2P", markevery=0.2)
    f2p_axis_rel[x_counter, y_counter].set_xlim([0, xs_cedar_static[-1]])
    f2p_axis_rel[x_counter, y_counter] \
        .set_ylim([0, max(ys_sead_static_relative[y] * 1.1 for y in range((int)(len(ys_sead_static_relative) / 2),
                                                                          len(ys_sead_static_relative)))])
    f2p_axis_abs[0, 0].legend(loc="upper left")
    f2p_axis_rel[0, 0].legend(loc="upper left")


def f2p(cnt_size, h_method):
    xs_f2p_static = []
    with open('res/F2P_n' + str(cnt_size) + '_h' + str(h_method) + '.res') as fp:
        lines = fp.readlines()
        for line in lines:
            xs_f2p_static.append(int(line.split("=", 1)[1]))
        # call to plot the current plot
    ys_f2p_static_abs = absolute_resolution(xs_f2p_static)
    ys_f2p_static_relative = relative_resolution(xs_f2p_static, ys_f2p_static_abs)
    return xs_f2p_static, ys_f2p_static_abs, ys_f2p_static_relative


# function to read the external files and make them lists
def res_files_to_G():
    f2p_figure_abs, f2p_axis_abs = plt.subplots(2, 2)
    f2p_figure_rel, f2p_axis_rel = plt.subplots(2, 2)
    # read the F2P information from external files and plot them in two different figures
    counter = 0
    for n in range(13, 15):
        for h in range(1, 3):
            f2p_VS_rest(f2p(n, h), n, counter, h - 1, f2p_axis_abs, f2p_axis_rel)
        counter += 1
    plt.tight_layout()
    plt.show()


# in order to check the convergence of the relative resolution in the cedar method
def toy_example(delta=0.5, max_val=200, cnt_size=8, exp_size=3, h_method=1):
    toy_fig_cedar, toy_axis_cedar = plt.subplots(2)
    xs_cedar_static, ys_cedar_static = cedar(delta, max_val)
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    toy_axis_cedar[0].plot(xs_cedar_static, ys_cedar_static)
    toy_axis_cedar[0].set_title("CEDAR Absolute Resolution")
    toy_axis_cedar[1].plot(xs_cedar_static, ys_cedar_relative)
    toy_axis_cedar[1].set_title("CEDAR Relative Resolution")

    toy_fig_sead, toy_axis_sead = plt.subplots(2)
    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size, exp_size)
    toy_axis_sead[0].plot(xs_sead_static, ys_sead_static_abs)
    toy_axis_sead[0].set_title("SEAD Absolute Resolution")
    toy_axis_sead[1].plot(xs_sead_static, ys_sead_static_relative)
    toy_axis_sead[1].set_title("SEAD Relative Resolution")

    toy_fig_f2p, toy_axis_f2p = plt.subplots(2)
    xs_f2p_static, ys_f2p_static_abs, ys_f2p_static_relative = f2p(cnt_size, h_method)
    toy_axis_f2p[0].plot(xs_f2p_static, ys_f2p_static_abs)
    toy_axis_f2p[0].set_title("F2P Absolute Resolution")
    toy_axis_f2p[1].plot(xs_f2p_static, ys_f2p_static_relative)
    toy_axis_f2p[1].set_title("F2P Relative Resolution")
    plt.tight_layout()
    plt.show()
