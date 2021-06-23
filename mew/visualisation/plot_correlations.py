#!/usr/bin/env python

from sys import argv

from mew.visualisation.plot import plot_correlation_vs_window

if __name__ == "__main__":
    in_dir = argv[1]
    out_dir = argv[2]
    window_size = int(argv[3])
    plot_correlation_vs_window(in_dir, out_dir, window_size)