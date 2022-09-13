import sys,json,time

from Bio import SeqIO
import random
from copy import deepcopy

def data_struct():
    return deepcopy({
        'generation':None,
        'seq_id':0,
        'is_alive':True,
        'parent_id':None,
        'path':'root',
        'mutations':{}
    })

def read_fasta_dict(fasta_file):
    """

    :param fasta_file: [str] Path to fasta file to read
    :return: [dict] of sequences indexed by sequence id
    """
    seqs = dict()
    with open(fasta_file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            seqs[str(record.id)] = str(record.seq).upper()
    handle.close()
    return seqs

def is_success(probability):
    if random.random() <= probability:
        return True
    else:
        return False

def get_available_positions(sites):
    pos = []
    for idx, bases in enumerate(sites):
        if sum(bases.values()) > 1:
            continue
        pos.append(idx)
    return pos


def simulate(seq,num_mutations,num_generations,num_children=2):
    seq_len = len(seq)
    seq_stack = [data_struct()]
    base_range = range(0,seq_len)
    bases = ['A','T','C','G']
    sites = []
    for i in base_range:
        sites.append( {} )
        for b in bases:
            sites[i][b] = 0
        sites[i][seq[i]] = 1

    available_sites = get_available_positions(sites)
    children = range(0,num_children)
    sample_id = 1
    events = range(0,num_mutations)
    for i in range(0,num_generations):
        if len(available_sites) == 0:
            print("Sequence has no more sites to mutate")
            break
        print("generation:{}, num seqs:{}".format(i,sample_id))
        for k in range(0,len(seq_stack)):
            parent_seq = seq_stack[k]
            if not parent_seq['is_alive']:
                continue
            parent_seq['is_alive'] = False

            for j in children:
                mutations = deepcopy(parent_seq['mutations'])
                for l in events:
                    pos = random.choice(available_sites)
                    muts = sites[pos]
                    for b in bases:
                        if muts[b] == 0:
                            break
                sites[pos][b] = 1
                available_sites = list(set(available_sites) - set([pos]))
                mutations[pos] = b
                child = data_struct()
                child['path'] = "{}.{}".format(parent_seq['path'],k)
                child['generation'] = i
                child['parent_id'] = k
                child['seq_id'] = sample_id
                child['mutations'] = mutations
                seq_stack.append(child)
                sample_id+=1
        available_sites = get_available_positions(sites)
    return seq_stack


def create_seqs(sim_samples,id,seq,outfile):
    fh = open(outfile,'w')
    fh.write(">{}\n{}\n".format(id,seq))
    for i in range(0,len(sim_samples)):
        data = sim_samples[i]
        if data['is_alive'] == False:
            continue
        id = "{}.{}".format(data['path'],data['seq_id'])
        res = list(seq)
        for pos in data['mutations']:
            res[pos] = data['mutations'][pos]

        fh.write(">{}\n{}\n".format(id,''.join(res)))
    fh.close()


def main():
    input_fasta = '/Users/jrobertson/Desktop/2022-09-Simulation/thrA.fas'
    output_fasta = '/Users/jrobertson/Desktop/2022-09-Simulation/thrA.simulation.fas'
    num_generations = 10
    inclusion_prob = 0.5
    num_mutations = 1


    fasta_seq = read_fasta_dict(input_fasta)
    ref_id = list(fasta_seq.keys())[0]
    ref_seq = fasta_seq[ref_id]


    sim_samples = simulate(ref_seq,num_mutations,num_generations)
    for i in range(0, len(sim_samples)):
        state = False
        num = random.uniform(0, 1)
        if num > inclusion_prob:
            state = True
        sim_samples[i]['is_alive'] = state

    create_seqs(sim_samples, ref_id,ref_seq, output_fasta)



if __name__ == '__main__':
    main()