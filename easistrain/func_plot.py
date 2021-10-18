# -*- coding: utf-8 -*-
"""
Created on Thu May 20 20:08:49 2021

@author: slim
"""

import matplotlib.pyplot as plt
import numpy as np


def showplot(x, y, xlabel, ylabel, pt, legend, name, title):
    plt.figure(num=name, figsize=(10, 8))
    for i in range(np.shape(x)[0]):
        plt.plot(x[i], y[i], pt, label=legend[i])
    plt.xlabel(xlabel, family="sans-serif", fontsize=28)
    plt.ylabel(ylabel, family="sans-serif", fontsize=28)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(loc="best", fontsize=22)
    plt.grid()
    plt.title(title, fontsize=30)
    plt.savefig(name, dpi=200)
    plt.close()
    return


# def saveplot(x,y,path)
