#!/usr/bin/env python
import argparse
from mew.visualisation.plot import *

def make_parser():
    parser = argparse.ArgumentParser(description="Plot multiple feature importance plots into a single figure.")

    parser.add_argument('-fi', '--feature_importances_files', nargs='+', type=str, required=True, help="Directories to .txt files containing feature importances")
    parser.add_argument('-en', '--encodings', nargs='+', type=str, required=True, help="Featurisation used to generate feature importances files")
    parser.add_argument('-l', '--labels', nargs='+', type=str, required=True, help="Labels for feature importances files")
    parser.add_argument('-t', '--title', type=str, required=True, help="Figure title")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output file")

    return parser

if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    plot_multiple(args.feature_importances_files, args.encodings, args.labels, args.title, args.output)
