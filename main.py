#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Libraries and files
import matplotlib.pyplot as plt

# Global variables
DOWN_TKT = "ListDown.tkt"
RI_TKT = "ListRI.tkt"
UP_TKT = "ListUp.tkt"


def read_file(file_name):
    temp_lst = []
    with open(file_name) as f: # file opening
        content = f.read()
        lines = content.split('\n')
        for l in lines:
            temp_lst.append(float(l))
    return temp_lst


def draw_plot_graph(l1x_values, l1y_values, title):
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 1, 1)
    plt.plot(x_values, y_values, label='Significant', color="red")
    plt.grid()
    plt.title(title)
    plt.xlabel("RI")
    plt.ylabel("Number of G4")
    plt.savefig("plot_RI.png")
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
    #down_exon = read_file(DOWN_TKT)
    ri_exon = read_file(RI_TKT)
    #up_exon = read_file(UP_TKT)

    x_values = list(range(len(ri_exon)))
    y_values = ri_exon

    draw_plot_graph(x_values, y_values, title="RI")

    # values = []
    # for i in range(len(x_values)):
    #     count_at_position = y_values[i]
    #
    #     for count in range(count_at_position):
    #         values.append(x_values[i])
    #
    # draw_hist_graph(values, title="RI", bins=len(x_values) // 100)
