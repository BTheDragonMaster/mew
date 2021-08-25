#!/usr/bin/env python

import argparse
from mew.visualisation.plot import plot_multiple_correlations

def make_parser():
    parser = argparse.ArgumentParser(description="Plot multiple feature importance plots into a single figure.")

    parser.add_argument('-cv', '--crossval_dirs', nargs='+', type=str, required=True, help="Directories to window crossval analyses")
    parser.add_argument('-w', '--window_sizes', nargs='+', type=int, required=True, help="Window sizes used to generate crossval analyses")
    parser.add_argument('-l', '--labels', nargs='+', type=str, required=True, help="Labels for crossval analyses")
    parser.add_argument('-t', '--title', type=str, required=True, help="Figure title")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file")
    parser.add_argument('-co', '--correlation_type', type=str, default='pearson', help="Type of correlation to plot (pearson or spearman)")

    return parser

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    plot_multiple_correlations(args.crossval_dirs, args.window_sizes, args.labels, args.output, args.title, correlation_type=args.correlation_type)