import os
from sys import argv
from pprint import pprint

from mew.parsers import read_correlations_per_group_from_dir


def get_best_models(models_dir):
    model_to_windows = {}
    for windows_dir_name in os.listdir(models_dir):
        windows_dir = os.path.join(models_dir, windows_dir_name)

        if os.path.isdir(windows_dir) and 'windows' in windows_dir_name:
            model_info = windows_dir_name.split('_')
            model = None
            featurisation = None
            window_size = None
            parameters = None
            if len(model_info) == 3:
                model, featurisation, window_size = model_info
            elif len(model_info) == 4:
                model, parameters, featurisation, window_size = model_info

            if model:
                window_size = int(window_size.split('windows')[1])
                if parameters:
                    alpha = float(parameters.split('alpha')[1])
                    model_tuple = (model, alpha, featurisation, window_size)
                else:
                    model_tuple = (model, featurisation, window_size)

                print(model_tuple)

                window_to_correlation = read_correlations_per_group_from_dir(windows_dir)
                model_to_windows[model_tuple] = window_to_correlation

    pprint(model_to_windows.keys())



    model_parameters_to_best_window_and_correlation = {}
    for model_parameters, windows_to_correlation in model_to_windows.items():
        best_pearson = -1.0
        best_pearson_stdev = None
        best_window = None
        for window, correlations in windows_to_correlation.items():
            pearson = correlations['pearson']
            if pearson > best_pearson:
                best_pearson = pearson
                best_window = window
                best_pearson_stdev = correlations['pearson_stdev']

        model_parameters_to_best_window_and_correlation[model_parameters] = (best_window, best_pearson, best_pearson_stdev)

    return model_parameters_to_best_window_and_correlation

def get_best_model(model_to_top_model):
    best_pearson = -1.0
    best_model = None
    best_window = None
    for model, top_model in model_to_top_model.items():
        window, pearson, pearson_stdev = top_model
        if pearson > best_pearson:
            best_pearson = pearson
            best_model = model
            best_window = window

    return best_model, best_window, best_pearson

def sort_models(model_to_top_model):
    models = list(model_to_top_model.keys())
    models.sort(key=lambda x: model_to_top_model[x][1], reverse=True)

    return models


def write_sorted_models(model_to_top_model, sorted_models, out_file):
    with open(out_file, 'w') as out:
        out.write('model\tfeaturisation\talpha\twindow size\tbest window\tpearson\tpearson_stdev\n')
        for model in sorted_models:

            if len(model) == 3:
                out.write(f'{model[0]}\t{model[1]}\tn/a\t{model[2]}\t')
            elif len(model) == 4:
                out.write(f'{model[0]}\t{model[2]}\t{model[1]}\t{model[3]}\t')

            window, pearson, pearson_stdev = model_to_top_model[model]
            out.write(f'{window}\t')
            out.write(format(pearson, '.4f'))
            out.write('\t')
            out.write(format(pearson_stdev, '.4f'))
            out.write('\n')


if __name__ == "__main__":
    top_models = get_best_models(argv[1])
    top_model = get_best_model(top_models)
    print(top_model)
    sorted_models = sort_models(top_models)
    write_sorted_models(top_models, sorted_models, argv[2])








