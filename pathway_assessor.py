import pandas as pd


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


