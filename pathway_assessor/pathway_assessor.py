import csv
import pathlib
import pickle
import os

from collections.abc import Iterable

import numpy as np
import pandas as pd
import scipy.stats as stats


def processed_expression_table(df):
    df.index.name = 'genes'
    return df.groupby('genes').mean()


def expression_ranks(expression_table_df, ascending, rank_method='max'):
    return expression_table_df.rank(method=rank_method, ascending=ascending)


# bg_genes: df of samples with background gene count
def bg_genes(expression_ranks_df):
    return expression_ranks_df.count()


def pathway_ranks(pathway_genes, expression_ranks_df, rank_method):
    return expression_ranks_df.reindex(pathway_genes).rank(method=rank_method).dropna(how='all')


def effective_pathway(pathway_ranks_df):
    return pathway_ranks_df.max()


def b(expression_ranks_df, pathway_ranks_df):
    return expression_ranks_df.subtract(pathway_ranks_df).dropna(how='all')


def c(effective_pathway_series, pathway_ranks_df):
    return effective_pathway_series - pathway_ranks_df


def d(bg_series, pathway_ranks_df, b_df, c_df):
    return bg_series - pathway_ranks_df - b_df - c_df


def sample_2x2(pathway_ranks_dict, b_dict, c_dict, d_dict):
    final_dict = {
        sample: {
            gene: [
                [val, b_dict[sample][gene]],
                [c_dict[sample][gene], d_dict[sample][gene]]
            ]
            for (gene, val) in genes.items()
        }
        for (sample, genes) in pathway_ranks_dict.items()
    }
    return pd.DataFrame(final_dict)


def clean_fisher_exact(table):
    try:
        if np.isnan(table).any():
            return np.nan
        else:
            return stats.fisher_exact(table, alternative='greater')[1]
    except ValueError:
        print(table)


def p_values(sample_2x2_df):
    return sample_2x2_df.apply(np.vectorize(clean_fisher_exact))


def neg_log(table):
    return -np.log(table)


def harmonic_average(iterable):
    if 0 in iterable:
        return 0

    reciprocal_iterable = [1/el for el in iterable if ~np.isnan(el)]
    denom = sum(reciprocal_iterable)

    if denom == 0:
        return np.nan
    else:
        return len(reciprocal_iterable) / denom


def geometric_average(iterable):
    try:
        clean_iterable = [el for el in iterable if ~np.isnan(el)]
        if not len(clean_iterable):
            return np.nan
        return np.exp(np.sum(np.log(clean_iterable)) / len(clean_iterable))
    except ZeroDivisionError:
        return 0


def user_pw_metadata_f(pw_data, output_dir_path):
    output_loc = '{}/user_pathways.tsv'.format(output_dir_path)
    pd.DataFrame.from_dict(pw_data, orient='index').to_csv(output_loc, sep='\t')


def pw_metadata_f(pw_db_choice, output_dir_path):
    output_loc = '{}/{}.tsv'.format(output_dir_path, pw_db_choice)
    pw_data = pickle.load(open('databases/metadata/{}.pkl'.format(pw_db_choice), 'rb'))
    pw_data.to_csv(output_loc, sep='\t')


def output_dir(output_dir_path):
    pathlib.Path(output_dir_path).mkdir(parents=True, exist_ok=True)


def user_pathways(f):
    pathway_db = {}
    pw_data = {}
    with open(f, 'r') as csv_in:
        reader = csv.reader(csv_in)
        for row in reader:
            pw = row[0]
            db = row[1]
            genes = set(row[2:])
            pathway_db[pw] = genes
            pw_data[pw] = {
                'db': db,
                'count': len(genes)
            }
    return pathway_db, pw_data


def validate_db_name(db_name):
    available_dbs = ['kegg', 'hallmark', 'reactome', 'hmdb_smpdb']
    if db_name.lower() not in available_dbs:
        raise ValueError(
            "{} not recognized. Available dbs: {}".format(db_name, ",".join(available_dbs))
        )
    return True


def db_pathways_dict(db_name):
    validate_db_name(db_name)
    db_parent = os.path.dirname(os.path.abspath(__file__))
    with open('{}/databases/{}.pkl'.format(db_parent, db_name.lower()), 'rb') as f:
        pathways = pickle.load(f)
    return pathways


def validate_pathways(pw_dict):
    if not isinstance(pw_dict, dict):
        raise TypeError("Pathways should be a dictionary of lists or sets")
    if any(not isinstance(gene_list, Iterable) for gene_list in pw_dict.values()):
        raise TypeError("Pathways should be a dictionary of lists or sets")

    return True


def all(
        expression_table,
        pathways=None,
        db='kegg',
        geometric=True,
        min_p_val=True,
        ascending=True,
        rank_method='max'
):

    if not pathways:
        pathways = db_pathways_dict(db)
    else:
        validate_pathways(pathways)

    harmonic_averages = [None] * len(pathways)
    geometric_averages = []
    min_p_vals = []

    if geometric:
        geometric_averages = [None] * len(pathways)
    if min_p_val:
        min_p_vals = [None] * len(pathways)

    expression_table_df = processed_expression_table(expression_table)
    expression_ranks_df = expression_ranks(expression_table_df, ascending=ascending, rank_method=rank_method)
    bg_genes_df = bg_genes(expression_ranks_df)

    sample_order = expression_table_df.columns

    # perform analysis for each pathway
    for i, pathway in enumerate(pathways):
        print('starting: {}'.format(pathway))
        pathway_ranks_df = pathway_ranks(pathways[pathway], expression_ranks_df, rank_method=rank_method)
        effective_pathway_df = effective_pathway(pathway_ranks_df)
        b_df = b(expression_ranks_df, pathway_ranks_df)
        c_df = c(effective_pathway_df, pathway_ranks_df)
        d_df = d(bg_genes_df, pathway_ranks_df, b_df, c_df)

        sample_2x2_df = sample_2x2(
            pathway_ranks_df.to_dict(),
            b_df.to_dict(),
            c_df.to_dict(),
            d_df.to_dict()
        )
        p_values_df = p_values(sample_2x2_df)

        # Harmonic averaging is default
        harmonic_averages_series = neg_log(p_values_df.apply(harmonic_average).loc[sample_order])
        harmonic_averages_series.name = pathway
        harmonic_averages[i] = harmonic_averages_series

        if geometric:
            geometric_averages_series = neg_log(p_values_df.apply(geometric_average).loc[sample_order])
            geometric_averages_series.name = pathway
            geometric_averages[i] = geometric_averages_series
        if min_p_val:
            min_p_vals_series = neg_log(p_values_df.min().loc[sample_order])
            min_p_vals_series.name = pathway
            min_p_vals[i] = min_p_vals_series
        print('finished: {}'.format(pathway))

    harmonic_averages_df = pd.concat(harmonic_averages, axis=1).T

    if geometric:
        geometric_averages_df = pd.concat(geometric_averages, axis=1).T
    else:
        geometric_averages_df = None

    if min_p_val:
        min_p_vals_df = pd.concat(min_p_vals, axis=1).T
    else:
        min_p_vals_df = None

    return {
               'harmonic': harmonic_averages_df,
               'geometric': geometric_averages_df,
                'min_p_val': min_p_vals_df
    }


def pa_stats(
        expression_table,
        mode='harmonic',
        pathways=None,
        db='kegg',
        ascending=True,
        rank_method='max'
):
    if not pathways:
        pathways = db_pathways_dict(db)
    else:
        validate_pathways(pathways)

    averages = [None] * len(pathways)

    expression_table_df = processed_expression_table(expression_table)
    expression_ranks_df = expression_ranks(expression_table_df, ascending=ascending, rank_method=rank_method)
    bg_genes_df = bg_genes(expression_ranks_df)

    sample_order = expression_table_df.columns

    # perform analysis for each pathway
    for i, pathway in enumerate(pathways):
        pathway_ranks_df = pathway_ranks(pathways[pathway], expression_ranks_df, rank_method=rank_method)
        effective_pathway_df = effective_pathway(pathway_ranks_df)
        b_df = b(expression_ranks_df, pathway_ranks_df)
        c_df = c(effective_pathway_df, pathway_ranks_df)
        d_df = d(bg_genes_df, pathway_ranks_df, b_df, c_df)

        sample_2x2_df = sample_2x2(
            pathway_ranks_df.to_dict(),
            b_df.to_dict(),
            c_df.to_dict(),
            d_df.to_dict()
        )
        p_values_df = p_values(sample_2x2_df)

        if mode == 'geometric':
            averages_series = neg_log(p_values_df.apply(geometric_average).loc[sample_order])
        if mode == 'harmonic':
            averages_series = neg_log(p_values_df.apply(harmonic_average).loc[sample_order])
        if mode == 'min':
            averages_series = neg_log(p_values_df.min().loc[sample_order])

        averages_series.name = pathway
        averages[i] = averages_series

    averages_df = pd.concat(averages, axis=1).T

    return averages_df


# Try doing this with a decorator
def harmonic(
        expression_table,
        pathways=None,
        db='kegg',
        ascending=True,
        rank_method='max'
):
    return pa_stats(expression_table, 'harmonic', pathways, db, ascending, rank_method)


def geometric(
        expression_table,
        pathways=None,
        db='kegg',
        ascending=True,
        rank_method='max'
):
    return pa_stats(expression_table, 'geometric', pathways, db, ascending, rank_method)


def min_p_val(
        expression_table,
        pathways=None,
        db='kegg',
        ascending=True,
        rank_method='max'
):
    return pa_stats(expression_table, 'min', pathways, db, ascending, rank_method)
