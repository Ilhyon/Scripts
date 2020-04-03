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


def draw_plot_graph(x_values0, y_values0, x_values1, y_values1, v):
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 1, 1)
    plt.plot(x_values0, y_values0, label='Non-Significant', color="black")
    plt.plot(x_values1, y_values1, label='Significant', color="red")
    plt.grid()
    plt.title('Plot of the number of pG4 per site / tot site for the virus '+\
        v+'and the location RI')
    plt.xlabel('Position in RI')
    plt.ylabel("Percent of site with a G4")
    plt.savefig(v+"_RI.png")
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
    files = {'kunvRI' : ['List_RI_0_kunvRI.txt', 'List_RI_1_kunvRI.txt'],
            'sinvRI' : ['List_RI_0_sinvRI.txt', 'List_RI_1_sinvRI.txt'],
            'yvfRI' : ['List_RI_0_yvfRI.txt', 'List_RI_1_yvfRI.txt'],
            'zikvRI' : ['List_RI_0_zikvRI.txt', 'List_RI_1_zikvRI.txt']}
    for v in files:
        VLoc0 = read_file(files[v][0])
        VLoc1 = read_file(files[v][1])

        x_values0 = list(range(len(VLoc0)))
        y_values0 = VLoc0
        x_values1 = list(range(len(VLoc1)))
        y_values1 = VLoc1

        draw_plot_graph(x_values0, y_values0, x_values1, y_values1, v)

    # values = []
    # for i in range(len(x_values)):
    #     count_at_position = y_values[i]
    #
    #     for count in range(count_at_position):
    #         values.append(x_values[i])
    #
    # draw_hist_graph(values, title="RI", bins=len(x_values) // 100)
