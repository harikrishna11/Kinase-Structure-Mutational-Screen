
def get_whereclause(gene):

    gene_bits = gene.split("-")
    if len(gene_bits) == 1:
        where_clause = "gene = '%s'" % gene_bits[0]
    elif len(gene_bits) == 2:
        where_clause = "gene = '%s' and offset_canonical = '%s'" % tuple(gene_bits)
    else:
        raise Exception("Bad syntax: gene can contain at most one - character")

    return where_clause
