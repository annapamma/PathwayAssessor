import pandas as pd


def expression_table(f):
    df = pd.read_csv(f, sep='\t', header=0, index_col=0)
    df.index.name = 'genes'
    return df.groupby('genes').mean()


def expression_ranks(expression_table_df, ascending):
    return expression_table_df.rank(method='first', ascending=ascending)
