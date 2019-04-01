import csv
import pathlib
import pickle

import numpy as np
import pandas as pd
import scipy.stats as stats


def expression_table(f):
    df = pd.read_csv(f, sep='\t', header=0, index_col=0)
    df.index.name = 'genes'
    return df.groupby('genes').mean()


def expression_ranks(expression_table_df, ascending):
    return expression_table_df.rank(method='first', ascending=ascending)


# bg_genes: df of samples with background gene count
def bg_genes(expression_ranks_df):
    return expression_ranks_df.count()


def pathway_ranks(pathway_genes, expression_ranks_df):
    return expression_ranks_df.reindex(pathway_genes).rank().dropna(how='all')


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
    try:
        clean_iterable = [el for el in iterable if ~np.isnan(el)]
        return len(clean_iterable) / sum([1/el for el in clean_iterable])
    except ZeroDivisionError:
        return 0


def geometric_average(iterable):
    try:
        clean_iterable = [el for el in iterable if ~np.isnan(el)]
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
