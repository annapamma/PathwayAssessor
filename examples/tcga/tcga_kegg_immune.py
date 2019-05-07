import csv
import os

import pandas as pd

import pathway_assessor as pa


def pathways_dict(f):
    with open(f, 'r') as csv_in:
        csv_in = csv.reader(csv_in, delimiter='\t')
        return {rows[0]: rows[1:] for rows in csv_in}


def write_table(score_type):
    if ascending:
        direction = 'suppression'
    else:
        direction = 'activation'
    output_f = '{}/scores/{}_{}_{}_{}'.format(examples_dir, tumor, pw_name, direction, score_type)
    scores[score_type].to_csv(output_f, sep='\t')
    print('Finished: {}'.format(output_f))


if __name__ == '__main__':
    pa_home = os.getenv('PATHWAY_ASSESSOR_HOME')
    tumor = 'blca_normal_slim'
    pw_name = 'kegg_immune'
    db_dir = '{}/pathway_assessor/databases/as_tables'.format(pa_home)
    ascending = True

    examples_dir = os.path.dirname(os.path.abspath(__file__))
    expression_table_f = '{}/{}'.format(examples_dir, tumor)
    pathways_f = '{}/{}.tsv'.format(db_dir, pw_name)

    expression_table = pd.read_csv(expression_table_f, sep='\t', index_col=0)
    pathways = pathways_dict(pathways_f)

    scores = pa.all(expression_table=expression_table, ascending=ascending, pathways=pathways)

    for method in scores:
        write_table(method)
