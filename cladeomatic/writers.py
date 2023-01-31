from cladeomatic.version import __version__


def write_clade_snp_report(clade_data, out_file):
    #TO DO Write report
    return

def print_params(params,outfile):
    fh = open(outfile,'w')
    p = vars(params)
    fh.write("version\t{}\n".format(__version__))
    for name in p:
        value = p[name]
        fh.write("{}\t{}\n".format(name,value))
    fh.close()
    return

def write_kmers(kmers, out_file):
    '''

    Parameters
    ----------
    kmers
    out_file

    Returns
    -------

    '''
    header = "index\tkseq\ttarget_pos\ttarget_base\taln_start\taln_end\n"
    fh = open(out_file, 'w')
    fh.write(header)
    i = 0
    for pos in kmers:
        for base in kmers[pos]:
            if len(kmers[pos][base]) == 0:
                continue
            for kseq in kmers[pos][base]:
                out_string = [i, kseq, pos, base, kmers[pos][base][kseq]['aln_start'], kmers[pos][base][kseq]['aln_end']
                              ]
                i += 1
                fh.write("{}\n".format("\t".join([str(x) for x in out_string])))
    fh.close()

def write_node_report(clade_data, outfile):
    '''
    Writes the Node information data to a tsv file
    :param clade_data: dict Node information structure
    :param outfile: str
    :return:
    '''
    fh = open(outfile, 'w')

    header = ['clade_id', 'num_members','pos', 'base', 'ave_dist', 'is_valid','is_selected', 'closest_clade_id', 'closest_sample_id',
              'closest_sample_dist', 'field_name', 'fisher_oddsr', 'fisher_p']
    fh.write("{}\n".format("\t".join([str(x) for x in header])))

    for clade_id in clade_data:
        num_bases = len(clade_data[clade_id]['bases'])
        for i in range(0, num_bases):
            row = [
                clade_id,
                clade_data[clade_id]['num_members'],
                clade_data[clade_id]['pos'][i]+1,
                clade_data[clade_id]['bases'][i],
                clade_data[clade_id]['ave_within_clade_dist'],
                clade_data[clade_id]['is_valid'],
                clade_data[clade_id]['is_selected'],
                clade_data[clade_id]['closest_clade_id'],
                clade_data[clade_id]['closest_sample_id'],
                clade_data[clade_id]['closest_sample_dist'],
                '',
                '',
                '',
            ]
            if len(clade_data[clade_id]['fisher']) == 0:
                fh.write("{}\n".format("\t".join([str(x) for x in row])))
            else:
                for field_name in clade_data[clade_id]['fisher']:
                    row[10] = field_name
                    row[11] = clade_data[clade_id]['fisher'][field_name]['oddsr']
                    row[12] = clade_data[clade_id]['fisher'][field_name]['p']
                    fh.write("{}\n".format("\t".join([str(x) for x in row])))

    fh.close()

    return



def write_snp_report(snp_data, outfile):
    '''
    Accepts snp_data dict data structure and writes snp details
    :param snp_data: dict() {[chrom_id] : {[position]: {[base]: :dict()}}}
    :param outfile: str
    :return:
    '''
    fh = open(outfile, 'w')
    fh.write(
        "chrom\tpos\tbase\tclade_id\tis_canonical\tnum_positive_clade_members\tnum_total_clade_members\tis_ref\toddsr\tp_value\n")
    for chrom in snp_data:
        for pos in snp_data[chrom]:
            for base in snp_data[chrom][pos]:
                row = [chrom, pos + 1, base,
                       snp_data[chrom][pos][base]['clade_id'],
                       snp_data[chrom][pos][base]['is_canonical'],
                       snp_data[chrom][pos][base]['num_clade_members'],
                       snp_data[chrom][pos][base]['num_members'],
                       snp_data[chrom][pos][base]['is_ref'],
                       snp_data[chrom][pos][base]['oddsr'],
                       snp_data[chrom][pos][base]['p_value']]
                fh.write("{}\n".format("\t".join([str(x) for x in row])))
    fh.close()
    return

def write_genotypes(genotypes, outfile):
    '''
    Accepts a list of sample genotype heirarchies and writes them to a file
    :param genotypes: dict of leaf names with heirarchal list of node id's to represent tree structure
    :param outfile: str
    :return:
    '''
    fh = open(outfile, 'w')
    fh.write("sample_id\tgenotype\n")
    for sample_id in genotypes:
        genotype = genotypes[sample_id]
        fh.write("{}\t{}\n".format(sample_id, genotype))
    fh.close()
    return

def write_scheme(header,scheme,outfile):
    fh = open(outfile, 'w')
    fh.write("{}\n".format("\t".join(header)))
    for i in range(0, len(scheme)):
        row = "\t".join([str(x) for x in list(scheme[i].values())])
        fh.write("{}\n".format(row))
    fh.close()