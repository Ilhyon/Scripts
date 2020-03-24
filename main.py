#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Libraries and files
import matplotlib.pyplot as plt

def read_file(file_name):
    temp_lst = []
    with open(file_name) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            temp_lst.append(float(l))
    return temp_lst


def draw_plot_graph(x_values0, y_values0, x_values1, y_values1, v, loc):
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 1, 1)
    plt.plot(x_values0, y_values0, label='Non-Significant', color="black")
    plt.plot(x_values1, y_values1, label='Significant', color="red")
    plt.grid()
    plt.title('Plot of the number of pG4 per site / tot site for the virus '+\
        v+'and the location '+loc)
    plt.xlabel('Position in ' + loc)
    plt.ylabel("Percent of site with a G4")
    plt.savefig(v+'_'+loc+".svg", format="svg")
    plt.close()


def draw_hist_graph(values, title, bins):
    plt.figure(figsize=(10, 5))
    ax = plt.hist(values, bins=bins, color="cyan", edgecolor="black")
    ax.yaxis.grid()
    plt.title(title)
    plt.xlabel("Localisation")
    plt.ylabel("Number of G4")
    plt.savefig("hist_RI.png")
    plt.close()


if __name__ == '__main__':
    files = {'kunvRI' : {'Up' : ['List_Up_0_kunvRI.txt', 'List_Up_1_kunvRI.txt'],
                        'Down' : ['List_Down_0_kunvRI.txt', 'List_Down_1_kunvRI.txt'],
                        'RI' : ['List_RI_0_kunvRI.txt', 'List_RI_1_kunvRI.txt']},
            'sinvRI' : {'Up' : ['List_Up_0_sinvRI.txt', 'List_Up_1_sinvRI.txt'],
                        'Down' : ['List_Down_0_sinvRI.txt', 'List_Down_1_sinvRI.txt'],
                        'RI' : ['List_RI_0_sinvRI.txt', 'List_RI_1_sinvRI.txt']},
            'yvfRI' : {'Up' : ['List_Up_0_yvfRI.txt', 'List_Up_1_yvfRI.txt'],
                        'Down' : ['List_Down_0_yvfRI.txt', 'List_Down_1_yvfRI.txt'],
                        'RI' : ['List_RI_0_yvfRI.txt', 'List_RI_1_yvfRI.txt']},
            'zikvRI' : {'Up' : ['List_Up_0_zikvRI.txt', 'List_Up_1_zikvRI.txt'],
                        'Down' : ['List_Down_0_zikvRI.txt', 'List_Down_1_zikvRI.txt'],
                        'RI' : ['List_RI_0_zikvRI.txt', 'List_RI_1_zikvRI.txt']}}
    for v in files:
        for loc in files[v]:
            VLoc0 = read_file(files[v][loc][0])
            VLoc1 = read_file(files[v][loc][1])

            x_values0 = list(range(len(VLoc0)))
            y_values0 = VLoc0
            x_values1 = list(range(len(VLoc1)))
            y_values1 = VLoc1

            draw_plot_graph(x_values0, y_values0, x_values1, y_values1, v, loc)

    # values = []
    # for i in range(len(x_values)):
    #     count_at_position = y_values[i]
    #
    #     for count in range(count_at_position):
    #         values.append(x_values[i])
    #
    # draw_hist_graph(values, title="RI", bins=len(x_values) // 100)
