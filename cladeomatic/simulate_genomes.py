import sys
from Bio import SeqIO
import random
from copy import deepcopy
from argparse import (ArgumentParser, ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter)

def parse_args():
    class CustomFormatter(ArgumentDefaultsHelpFormatter, RawDescriptionHelpFormatter):
        pass

    parser = ArgumentParser(
        description="Generate generational mutated sequences",
        formatter_class=CustomFormatter)
    parser.add_argument('--in_fasta', type=str, required=False,
                        help='Fasta file to mutate',default=None)
    parser.add_argument('--num_generations', type=int, required=False,
                        help='Number of generations to simulate',default=4)
    parser.add_argument('--num_children', type=int, required=False,
                        help='Number of children per sequence per generation',default=5)
    parser.add_argument('--min_num_mutations', type=int, required=False,
                        help='Number of mutations per sequence per generation',default=1)
    parser.add_argument('--max_num_mutations', type=int, required=False,
                        help='Number of mutations per sequence per generation',default=5)
    parser.add_argument('--inclusion_prob', type=float, required=False,
                        help='Probability that a given sequence will be included in the final output',default=0.75)
    parser.add_argument('--length', type=int, required=False,
                        help='Length of random sequence',default=None)
    parser.add_argument('--outfile', type=str, required=True, help='Output file to put results')


    return parser.parse_args()

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


def simulate(seq,min_num_mutations,max_num_mutations,num_generations,num_children=2):
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

    for i in range(0,num_generations+1):
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
                if min_num_mutations != max_num_mutations:
                    num_mutations = random.randrange(min_num_mutations,max_num_mutations+1)
                else:
                    num_mutations = min_num_mutations

                events = range(0, num_mutations)
                if len(available_sites) == 0:
                    break
                for l in events:
                    pos = random.choice(available_sites)
                    muts = sites[pos]
                    random.shuffle(bases)
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


def create_seqs(sim_samples,seq,outfile):
    fh = open(outfile,'w')
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

def generate_random_sequence(length,a=0.25,t=0.25,c=0.25,g=0.25):
    seq = []
    total = a+t+c+g
    if total != 1:
        print("Error base composisition must sum to 1")
        sys.exit()

    for i in range(0,length):
        num = random.uniform(0, 1)
        if num >= 0 and num < a:
            b = 'A'
        elif num >= a and num < a + t:
            b = 'T'
        elif num >= a+t and num < a + t + c:
            b = 'C'
        else:
            b = 'G'
        seq.append(b)

    return ''.join(seq)

def main():
    cmd_args = parse_args()
    input_fasta = cmd_args.in_fasta
    output_fasta = cmd_args.outfile
    length = cmd_args.length
    num_children = cmd_args.num_children
    num_generations = cmd_args.num_generations
    inclusion_prob = 1 - cmd_args.inclusion_prob
    min_num_mutations = cmd_args.min_num_mutations
    max_num_mutations = cmd_args.max_num_mutations

    if input_fasta is None:
        if length is None:
            count_seqs = 1
            for i in range(0,num_generations):
                tmp = 0
                for j in range(0,count_seqs):
                    for k in range(0,num_children):
                        tmp+=1
                count_seqs = tmp
            length = count_seqs * max_num_mutations
        ref_seq = generate_random_sequence(length*10)
    else:
        fasta_seq = read_fasta_dict(input_fasta)
        ref_id = list(fasta_seq.keys())[0]
        ref_seq = fasta_seq[ref_id]

    sim_samples = simulate(ref_seq,min_num_mutations,max_num_mutations,num_generations,num_children)
    for i in range(0, len(sim_samples)):
        state = False
        num = random.uniform(0, 1)
        if num > inclusion_prob:
            state = True
        sim_samples[i]['is_alive'] = state

    create_seqs(sim_samples, ref_seq, output_fasta)



if __name__ == '__main__':
    main()