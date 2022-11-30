import math
import matplotlib
import matplotlib.pyplot as plt
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
CNT_SIZE = 8


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


def absolute_resolution(x_arr):
    y_arr = [1]
    for c in range(1, len(x_arr)):
        y_arr.append(x_arr[c] - x_arr[c - 1])
    return y_arr


def relative_resolution(x_arr, y_abs):
    y_arr = [1]
    for c in range(1, len(x_arr)):
        y_arr.append(y_abs[c] / x_arr[c])
    return y_arr


def dynamic_pre_calculate_stages(cnt_size):
    expansion_array = [(2 ** exp) for exp in range(0, cnt_size - 1)]
    stage = [0]
    for j in range(1, cnt_size - 1):
        stage.append(stage[j - 1] + (expansion_array[j] * 2 ** (cnt_size - 1 - j)))
    return stage


def dynamic_sead(cnt_size):
    """
    function that gets the number of bits,
    Returns lists of
    1. The counted numbers
    2. The absolute resolution
    3. The relative resolution
    in the Dynamic SEAD method
    """
    stage = dynamic_pre_calculate_stages(cnt_size)
    xs_dynamic = []
    for i in range((2 ** cnt_size) - 2):
        bin_cntr = np.binary_repr(i, cnt_size)
        # Counting the ones in the binary repr string, by looking fot the first 0 index from the left
        exp_size = bin_cntr.index('0')
        # Calculating the mantissa
        mantissa = int(bin_cntr[exp_size + 1:cnt_size], 2)
        xs_dynamic.append((mantissa * (2 ** exp_size)) + stage[exp_size])
    ys_dynamic_abs = absolute_resolution(xs_dynamic)
    ys_dynamic_relative = relative_resolution(xs_dynamic, ys_dynamic_abs)
    return xs_dynamic, ys_dynamic_abs, ys_dynamic_relative


def static_pre_calculate_stages(cnt_size, exp_size):
    expansion_array = [(2 ** exp) for exp in range(0, 2 ** exp_size)]
    expansion_array_sum = [1]
    for i in range(1, len(expansion_array)):
        expansion_array_sum.append(expansion_array_sum[i - 1] + expansion_array[i])
    stage = [((2 ** (cnt_size - exp_size)) * expansion_array_sum[j]) for j in range(0, len(expansion_array_sum) - 1)]
    stage.insert(0, 0)
    return stage, expansion_array


def static_sead(cnt_size, exp_size):
    """
     function that gets the number of bits and, the number of Exponent bits,
     Returns lists of
    1. the counted numbers
    2. the absolute resolution
    3. the relative resolution
    in the Static SEAD method
    """
    xs_static = []
    stage, expansion_array = static_pre_calculate_stages(cnt_size, exp_size)
    for e in range(2 ** exp_size):
        for m in range(2 ** (cnt_size - exp_size)):
            xs_static.append((m * 2 ** e) + stage[e])
    ys_static_abs = absolute_resolution(xs_static)
    ys_static_relative = relative_resolution(xs_static, ys_static_abs)
    return xs_static, ys_static_abs, ys_static_relative


# This is the CEDAR formula to calculate the diff given the delta and the sum of the previous diffs
calc_diff = lambda delta, sum_of_prev_diffs: (1 + 2 * delta ** 2 * sum_of_prev_diffs) / (1 - delta ** 2)


# function that computes the shared estimators and D by using the CEDAR's formula
def cedar(delta, max_val):
    shared_estimators = []
    different = []
    shared_estimators.append(0)
    i = 0
    while shared_estimators[i] < max_val:
        # using the cedar's formula
        different.append(calc_diff(delta=delta, sum_of_prev_diffs=shared_estimators[i]))
        i += 1
        shared_estimators.append(shared_estimators[i - 1] + different[i - 1])
    shared_estimators.pop()
    return shared_estimators, different


def ideal_exp_size(counted_num, cnt_size):
    """
    Calculate the min number of exponent bites, that you need in order to count up to 'counted_num',
    given the number of bites and the counted number
    raise value error if the counted num can't be reached
    """
    for ideal_exp in range(1, cnt_size):
        if ((2 ** (cnt_size - ideal_exp)) - 1) * 2 ** ((2 ** ideal_exp) - 1) >= counted_num:
            return ideal_exp
        # can't be reached
    return cnt_size - 1


def min_delta(init_delta, max_val, cnt_size, delta_jumps=DELTA_JUMPS):
    """
    calc the minimal delta that can be obtained for CEDAR, given the max_val and the cnt size
    This function assumes that the initial delta is suitable for the cnt_size and for the max value
    """
    for delta in np.arange(init_delta, 0, -delta_jumps):
        cur_val = 0
        for i in range(2 ** cnt_size):
            cur_val += calc_diff(delta=delta, sum_of_prev_diffs=cur_val)
        if cur_val < max_val:
            return round(delta + delta_jumps, 3)


def update_byte(val):
    update_graph1(CNT_SIZE, val)


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
    # Find fair conditions in order to compare between the different methods
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


def f2p_VS_rest(axis, method, cnt_size, x_counter, y_counter, f2p_axis_abs, f2p_axis_rel):
    xs_f2p, ys_f2p_abs, ys_f2p_relative = axis
    # for fare comparison with CEDAR we need to check the min delta
    min_del = min_delta(0.1, xs_f2p[-1], cnt_size, 0.005)
    # for fare comparison with CEDAR we need to check the min exp
    min_exp = ideal_exp_size(xs_f2p[-1], cnt_size)
    # compare to SEAD method and CEDAR to F2P
    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size, min_exp)
    xs_cedar_static, ys_cedar_static = cedar(min_del, xs_f2p[-1])
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    f2p_axis_abs[x_counter, y_counter].set_title("F" + str(method) + "P " +
        str(cnt_size) + " bits, " + str(round(min_del, 3)) + " Delta, " + str(
        min_exp) + " exp, " + "max " + str(
        xs_f2p[-1]))
    f2p_axis_abs[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_static, "-b", label="Static cedar")
    f2p_axis_abs[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_abs, "-r", label="Static SEAD")
    f2p_axis_abs[x_counter, y_counter].plot(xs_f2p, ys_f2p_abs, "-gD", label="F" + str(method) + "P", markevery=MARKERS_GAP)
    f2p_axis_abs[x_counter, y_counter].set_xlim([16, xs_cedar_static[-1]])
    f2p_axis_abs[x_counter, y_counter].set_ylim([0, ys_f2p_abs[-1] * 5])
    f2p_axis_abs[x_counter, y_counter].set_xlabel('Counted Value')
    f2p_axis_abs[x_counter, y_counter].set_ylabel('Absolute Resolution')
    plt.tight_layout()

    f2p_axis_rel[x_counter, y_counter].set_title("F" + str(method) + "P " +
        str(cnt_size) + " bits, " + str(round(min_del, 3)) + " Delta, " + str(min_exp) + " exp, " + "max " + str(
            xs_f2p[-1]))
    f2p_axis_rel[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_relative, "-b", label="Static cedar")
    f2p_axis_rel[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_relative, "-r", label="Static SEAD")
    f2p_axis_rel[x_counter, y_counter].plot(xs_f2p, ys_f2p_relative, "-gD", label="F"+ str(method)+"P", markevery=MARKERS_GAP)
    f2p_axis_rel[x_counter, y_counter].set_xlim([128, xs_cedar_static[-1]])
    f2p_axis_rel[x_counter, y_counter].set_yscale("log")

    # f2p_axis_rel[x_counter, y_counter] \
    #     .set_ylim([0, max(ys_sead_static_relative[y] * 1.01 for y in range((int)(len(ys_sead_static_relative) / 2),
    #                                                                       len(ys_sead_static_relative)))])
    f2p_axis_rel[x_counter, y_counter].set_ylim([0,0.1])
    f2p_axis_rel[x_counter, y_counter].set_xlabel('Counted Value')
    f2p_axis_rel[x_counter, y_counter].set_ylabel('Relative Resolution')
    f2p_axis_abs[x_counter, y_counter].legend(loc='upper left')
    f2p_axis_rel[x_counter, y_counter].legend(loc='best')


def f2p(method, cnt_size, h_method):
    xs_f2p_static = []
    with open('res/F' + str(method) + 'P_n' + str(cnt_size) + '_h' + str(h_method) + '.res') as fp:
        lines = fp.readlines()
        for line in lines:
            xs_f2p_static.append(int(line.split("=", 1)[1]))
        # call to plot the current plot
    ys_f2p_static_abs = absolute_resolution(xs_f2p_static)
    ys_f2p_static_relative = relative_resolution(xs_f2p_static, ys_f2p_static_abs)
    return xs_f2p_static, ys_f2p_static_abs, ys_f2p_static_relative


# Function to read the external files and make them lists
def read_f2p_files_to_4graphs(method: list, cnt_size: list, h_method: list):
    f2p_figure_abs, f2p_axis_abs = plt.subplots(2, 2)
    f2p_figure_rel, f2p_axis_rel = plt.subplots(2, 2)

    # read the F2P information from external files and plot them in two different figures
    for i in range(0, 4):
        f2p_VS_rest(f2p(method[i], cnt_size[i], h_method[i]), method[i], cnt_size[i], i % 2, (int)(i / 2), f2p_axis_abs,
                    f2p_axis_rel)
    plt.tight_layout()
    plt.show()


def read_one_file_to_graph(method, cnt_size, h_method):
    f2p_figure_abs, f2p_axis_abs = plt.subplots(2, 2)
    f2p_figure_rel, f2p_axis_rel = plt.subplots(2, 2)
    f2p_VS_rest(f2p(method, cnt_size, h_method), cnt_size, 0, 0, f2p_axis_abs, f2p_axis_rel)
    plt.tight_layout()
    plt.show()


def generate_txt_file(cnt_size, max_val, method="cedar"):
    newfile = open(str(method)+"_up_to_"+str(max_val)+"_using_"+str(cnt_size)+"_bites"+".txt", "w")
    min_del = min_delta(0.1, max_val, cnt_size)
    newfile.write("The min Delta: "+str(min_del) + '\n')
    xs = cedar(min_del, max_val)
    for x in xs[0]:
        newfile.write(str(x)+'\n')
    newfile.close()


# In order to check the convergence of the relative resolution in the different methods
def toy_example(delta=0.5, max_val=200, cnt_size=4, exp_size=1, h_method=1):
    toy_fig_cedar, toy_axis_cedar = plt.subplots(2)
    xs_cedar_static, ys_cedar_static = cedar(delta, max_val)
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    toy_axis_cedar[0].plot(xs_cedar_static, ys_cedar_static, "-gD")
    toy_axis_cedar[0].set_title("CEDAR Absolute Resolution")
    toy_axis_cedar[1].plot(xs_cedar_static, ys_cedar_relative, "-gD")
    toy_axis_cedar[1].set_title("CEDAR Relative Resolution")

    toy_fig_sead, toy_axis_sead = plt.subplots(2)
    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size, exp_size)
    toy_axis_sead[0].plot(xs_sead_static, ys_sead_static_abs, "-gD")
    toy_axis_sead[0].set_title("SEAD Absolute Resolution")
    toy_axis_sead[1].plot(xs_sead_static, ys_sead_static_relative, "-gD")
    toy_axis_sead[1].set_title("SEAD Relative Resolution")

    toy_fig_f2p, toy_axis_f2p = plt.subplots(2)
    xs_f2p_static, ys_f2p_static_abs, ys_f2p_static_relative = f2p(CNT_SIZE, h_method)
    toy_axis_f2p[0].plot(xs_f2p_static, ys_f2p_static_abs, "-gD", markevery=MARKERS_GAP)
    toy_axis_f2p[0].set_title("F2P Absolute Resolution")
    toy_axis_f2p[1].plot(xs_f2p_static, ys_f2p_static_relative, "-gD", markevery=MARKERS_GAP)
    toy_axis_f2p[1].set_title("F2P Relative Resolution")
    plt.tight_layout()
    plt.show()
