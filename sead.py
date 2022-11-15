import math
import matplotlib
import matplotlib.pyplot as plt
import numpy
from matplotlib.widgets import Slider

MARKER_SIZE = 16
MARKER_SIZE_SMALL = 1
LINE_WIDTH = 3
LINE_WIDTH_SMALL = 1
FONT_SIZE = 20
FONT_SIZE_SMALL = 5
LEGEND_FONT_SIZE = 14
LEGEND_FONT_SIZE_SMALL = 5
UNIFORM_CHAIN_MIG_COST = 600


# Set the parameters of the plot (sizes of fonts, legend, ticks etc.).
# mfc='none' makes the markers empty.
def set_plt_params(size='large'):
    if size == 'large':
        matplotlib.rcParams.update({'font.size': FONT_SIZE,
                                    'legend.fontsize': LEGEND_FONT_SIZE,
                                    'xtick.labelsize': FONT_SIZE,
                                    'ytick.labelsize': FONT_SIZE,
                                    'axes.labelsize': FONT_SIZE,
                                    'axes.titlesize': FONT_SIZE, })
    else:
        matplotlib.rcParams.update({'font.size': FONT_SIZE_SMALL,
                                    'legend.fontsize': LEGEND_FONT_SIZE_SMALL,
                                    'xtick.labelsize': FONT_SIZE_SMALL,
                                    'ytick.labelsize': FONT_SIZE_SMALL,
                                    'axes.labelsize': FONT_SIZE_SMALL,
                                    'axes.titlesize': FONT_SIZE_SMALL, })


def pre_calculate_stages(cnt):
    stage = [1]
    for i in range(1, 2 ** cnt):
        stage.append(stage[i - 1] + (2 ** cnt))
    return stage


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


def from_binary(arr, start, end):
    k = 0
    count_sum = 0
    for i in range(start, end + 1):
        if arr[i]:
            count_sum += pow(2, k)
        k += 1
    return count_sum


def binary_to_val_static_sead(arr, exp_size):
    mantissa = from_binary(arr, exp_size, len(arr) - 1)
    exponent = from_binary(arr, 0, exp_size - 1)
    return mantissa * 2 ** exponent


def binary_to_val_dynamic_sead(arr, stage):
    exp_d = 0
    while arr[exp_d]:
        exp_d += 1
        if exp_d == 8:
            return 1
    mantissa = from_binary(arr, exp_d + 1, len(arr) - 1)
    return mantissa * 2 ** exp_d + stage[exp_d]


def dynamic_sead(cnt_size):
    bits_arr = numpy.zeros(cnt_size, bool)
    xs_dynamic = []
    finish = True
    i = 0
    # filing the dynamic array
    stage = pre_calculate_stages(cnt_size)
    while finish:
        xs_dynamic.append(binary_to_val_dynamic_sead(bits_arr, stage))
        while True:
            i += 1
            if not bits_arr[i - 1]:
                break
            elif i == cnt_size:
                finish = False
                break
        i -= 1
        bits_arr[i] = True
        while i > 0:
            i -= 1
            bits_arr[i] = False
    xs_dynamic = sort_and_erase_duplicates(xs_dynamic)
    ys_dynamic_abs = absolute_resolution(xs_dynamic)
    ys_dynamic_relative = relative_resolution(xs_dynamic, ys_dynamic_abs)
    return xs_dynamic, ys_dynamic_abs, ys_dynamic_relative


def static_sead(cnt_size, exp_size):
    xs_static = []
    for e in range(2 ** exp_size):
        for m in range(2 ** (cnt_size - exp_size)):
            xs_static.append(m * 2 ** e)
    xs_static = sort_and_erase_duplicates(xs_static)
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
    return ideal_exp_size_recursive(counted_num, 1, cnt_size)


def ideal_exp_size_recursive(counted_num, exp_l, cnt_size):
    if exp_l == cnt_size - 1:
        return exp_l
    elif counted_num > ((2 ** (cnt_size - exp_l)) - 1) * 2 ** ((2 ** exp_l) - 1):
        return ideal_exp_size_recursive(counted_num, exp_l + 1, cnt_size)
    return exp_l


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
        return delta + 0.01
    return min_delta(delta - 0.01, max_val, cnt_size)


def sort_and_erase_duplicates(xs):
    xs.sort()
    temp_arr = []
    for x in xs:
        if x not in temp_arr:
            temp_arr.append(x)
    xs = temp_arr
    return xs


def update_byte(val):
    plot_graphs(8, val)


def delta_changed(val):
    update_cedar(val, max_value.val)


def max_value_changed(val):
    update_cedar(delta_slider.val, val)


def update_cedar(delta, max_val):
    plot_graph2(delta, max_val)


# plot the static and dynamic sead to one graph
def plot_graphs(cnt_size, exp_size):
    xs_static, ys_static_abs, ys_static_relative = static_sead(cnt_size, exp_size)
    xs_dynamic, ys_dynamic_abs, ys_dynamic_relative = dynamic_sead(cnt_size)
    axis[0, 0].clear()
    axis[0, 0].plot(xs_static, ys_static_abs)
    axis[0, 0].set_title('Static SEAD')
    axis[0, 0].set_xlabel('Real Value')
    axis[0, 0].set_ylabel('absolute resolution')
    axis[0, 0].set_xlim([0, xs_static[-1]])
    axis[0, 0].set_ylim(0)

    axis[1, 0].clear()
    axis[1, 0].plot(xs_static, ys_static_relative, 'C3')
    axis[1, 0].set_xlabel('Real Value')
    axis[1, 0].set_ylabel('relative resolution')
    axis[1, 0].set_ylim([0, 0.15])
    axis[1, 0].set_xlim([0, xs_static[-1]])

    axis[0, 1].clear()
    axis[0, 1].plot(xs_dynamic, ys_dynamic_abs)
    axis[0, 1].set_title('Dynamic SEAD')
    axis[0, 1].set_xlabel('Real Value')
    axis[0, 1].set_ylabel('absolute resolution')
    axis[0, 1].set_xlim([0, xs_dynamic[-1]])
    axis[0, 1].set_ylim(0)

    axis[1, 1].clear()
    axis[1, 1].plot(xs_dynamic, ys_dynamic_relative, 'C3')
    axis[1, 1].set_xlabel('Real Value')
    axis[1, 1].set_ylabel('relative resolution')
    axis[1, 1].set_ylim([0, 0.15])
    axis[1, 1].set_xlim([0, xs_dynamic[-1]])
    plt.tight_layout()
    plt.show()


# plot the CEDAR and the SEAD at the same graph
def plot_graph2(delta, max_val):
    xs_cedar_static, ys_cedar_static = cedar(delta, max_val)
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)

    # find fair conditions to compare between the different methods
    cnt_size_needed = math.ceil(math.log2(len(xs_cedar_static)))
    exp_size_needed = ideal_exp_size(xs_cedar_static[-1], cnt_size_needed)

    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size_needed, exp_size_needed)
    xs_sead_dynamic, ys_sead_dynamic_abs, ys_sead_dynamic_relative = dynamic_sead(cnt_size_needed)

    axs[0].clear()
    axs[0].plot(xs_cedar_static, ys_cedar_static, "-b", label="Static cedar")
    axs[0].plot(xs_sead_static, ys_sead_static_abs, "-r", label="Static SEAD")
    axs[0].plot(xs_sead_dynamic, ys_sead_dynamic_abs, "-g", label="Dynamic SEAD")
    axs[0].legend(loc="upper left")
    axs[0].set_title(str(exp_size_needed) + " Exp bites and " + str(cnt_size_needed) +
                     " bits need in order to count up to " + str(math.ceil(xs_cedar_static[-1])) + " ")
    axs[0].set_xlabel('Real Value')
    axs[0].set_ylabel('absolute resolution')
    axs[0].set_xlim([0, xs_cedar_static[-1]])
    axs[0].set_ylim(0, 120)

    axs[1].clear()
    axs[1].plot(xs_cedar_static, ys_cedar_relative, "-b", label="Static cedar")
    axs[1].plot(xs_sead_static, ys_sead_static_relative, "-r", label="Static SEAD")
    axs[1].plot(xs_sead_dynamic, ys_sead_dynamic_relative, "-g", label="Dynamic SEAD")
    axs[1].set_xlim([10, xs_cedar_static[-1]])
    axs[1].set_ylim(0, 0.18)
    axs[1].set_xlabel('Real Value')
    axs[1].set_ylabel('relative resolution')


def F2P_VS_rest(xs_f2p, cnt_size, x_counter, y_counter):
    # calculate the absolute resolution for f2p
    ys_f2p_abs = absolute_resolution(xs_f2p)
    # calculate the relative resolution for f2p
    ys_f2p_relative = relative_resolution(xs_f2p, ys_f2p_abs)
    min_del = min_delta(0.1, xs_f2p[-1], cnt_size)
    min_exp = ideal_exp_size(xs_f2p[-1], cnt_size)
    # compare to SEAD method and CEDAR
    xs_sead_static, ys_sead_static_abs, ys_sead_static_relative = static_sead(cnt_size, min_exp)
    xs_cedar_static, ys_cedar_static = cedar(min_del, xs_f2p[-1])
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)

    f2p_axis_abs[x_counter, y_counter].set_title(
        str(cnt_size) + " bits, " + str(round(min_del, 2)) + " Delta, " + str(min_exp) + " exp, " + "max " + str(
            xs_f2p[-1]))
    f2p_axis_abs[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_static, "-b", label="Static cedar")
    f2p_axis_abs[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_abs, "-r", label="Static SEAD")
    f2p_axis_abs[x_counter, y_counter].plot(xs_f2p, ys_f2p_abs, "-g", label="F2P")
    f2p_axis_abs[x_counter, y_counter].set_xlim([0, xs_cedar_static[-1]])
    f2p_axis_abs[x_counter, y_counter].set_ylim([0, ys_f2p_abs[-1] * 5])

    f2p_axis_rel[x_counter, y_counter].set_title(
        str(cnt_size) + " bits, " + str(round(min_del, 2)) + " Delta, " + str(min_exp) + " exp, " + "max " + str(
            xs_f2p[-1]))
    f2p_axis_rel[x_counter, y_counter].plot(xs_cedar_static, ys_cedar_relative, "-b", label="Static cedar")
    f2p_axis_rel[x_counter, y_counter].plot(xs_sead_static, ys_sead_static_relative, "-r", label="Static SEAD")
    f2p_axis_rel[x_counter, y_counter].plot(xs_f2p, ys_f2p_relative, "-g", label="F2P")
    f2p_axis_rel[x_counter, y_counter].set_xlim([0, xs_cedar_static[-1]])
    f2p_axis_rel[x_counter, y_counter].set_ylim([0,
                                                 max(ys_sead_static_relative[y] for y in
                                                     range((int)(len(ys_sead_static_relative) / 2),
                                                           len(ys_sead_static_relative)))])
    plt.tight_layout()


def res_files_to_G():
    x_counter = 0
    n_counter = 13
    while n_counter < 15:
        h_counter = 1
        while h_counter < 3:
            xs = []
            with open('res/F2P_n' + str(n_counter) + '_h' + str(h_counter) + '.res') as fp:
                lines = fp.readlines()
                for line in lines:
                    xs.append(int(line.split("=", 1)[1]))
                F2P_VS_rest(xs, n_counter, x_counter, h_counter - 1)
            h_counter += 1
        n_counter += 1
        x_counter += 1


# in order to check the convergence of the relative resolution in the cedar method
def toy_example(delta, max_val):
    xs_cedar_static, ys_cedar_static = cedar(delta, max_val)
    ys_cedar_relative = relative_resolution(xs_cedar_static, ys_cedar_static)
    toy_axis[1].plot(xs_cedar_static, ys_cedar_relative)
    toy_axis[1].set_title("relative")
    toy_axis[0].plot(xs_cedar_static, ys_cedar_static)
    toy_axis[0].set_title("abs")


# set_plt_params()
f2p_figure_abs, f2p_axis_abs = plt.subplots(2, 2)
f2p_figure_rel, f2p_axis_rel = plt.subplots(2, 2)

res_files_to_G()
f2p_axis_abs[0, 0].legend(loc="upper left")
f2p_axis_rel[0, 0].legend(loc="upper left")
figure, axis = plt.subplots(2, 2)
byte_slider = Slider(plt.axes([0.53, 0.95, 0.1, 0.04]), 'Exp', 2, 7, 3)
byte_slider.valstep = 1
byte_slider.on_changed(update_byte)

# second figure
fig, axs = plt.subplots(2)
figure.subplots_adjust(0.25, 0.25)
delta_slider = Slider(plt.axes([0.1, 0.95, 0.1, 0.04]), 'Delta', 0.01, 0.2, 0.1)
delta_slider.valstep = 0.01
max_value = Slider(plt.axes([0.32, 0.95, 0.1, 0.04]), 'Max Value', 255, 5000, 2048)
max_value.valstep = 20
plot_graph2(delta_slider.val, max_value.val)
delta_slider.on_changed(delta_changed)
max_value.on_changed(max_value_changed)
# toy example
toy_fig, toy_axis = plt.subplots(2)
toy_example(0.5, 200)
plot_graphs(8, 3)
