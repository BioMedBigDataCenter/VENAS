# __title__   = Find parsimony informative sites
# __created__ = 2020-09-09
# __author__  = Ruifang Cao
# __python__  = 3.+

from Bio import AlignIO, SeqIO
import os
from tqdm import tqdm
import numpy as np
import argparse
import sys

def find_ps_sites(alignment, seq, remove_fq, epis_base):
    rm_pis = []
    all_pis = []
    default_nc = {'a', 't', 'g', 'c', 'A', 'T', 'C', 'G'}
    for width in tqdm(range(0, len(alignment[0])), desc='First step: find ePIS'):
        uniq_nc = {}
        tmp_all = 0
        for j in range(0, len(alignment)):
            tmp_key = alignment[j].seq[width]
            if tmp_key not in default_nc:
                continue
            elif tmp_key in uniq_nc.keys():
                uniq_nc[tmp_key] = uniq_nc[tmp_key] + 1
                tmp_all += 1
            else:
                uniq_nc[tmp_key] = 1
                tmp_all += 1
        # The number of effective bases is greater than the parameter "base"
        if len(uniq_nc) > 1 and tmp_all > epis_base:
            cur = 0
            for key in uniq_nc.keys():
                if uniq_nc[key] > 1:
                    cur += 1

            if cur > 1 and seq[width] != '-':
                max_dic = max(uniq_nc, key=uniq_nc.get)
                pos_freq = 1 - (uniq_nc[max_dic] / tmp_all)
                all_pis.append(width)
                # remove_fq is a user-defined parameter for which frequencies need to be removed
                if pos_freq > remove_fq:
                    rm_pis.append(width)
    return all_pis, rm_pis


def pis_freq(align_filter, sites, seq, f_base):
    freqs = []
    default_nc = {'a', 't', 'g', 'c', 'A', 'T', 'C', 'G'}
    for width in tqdm(sites, desc='Second step: calculate ePIS frequency'):
        uniq_nc = {}
        tmp_all = 0
        tmp_freq = []
        for j in range(0, len(align_filter)):
            tmp_key = align_filter[j].seq[width]
            if tmp_key not in default_nc:
                continue
            elif tmp_key in uniq_nc.keys():
                uniq_nc[tmp_key] = uniq_nc[tmp_key] + 1
                tmp_all += 1
            else:
                uniq_nc[tmp_key] = 1
                tmp_all += 1
        # rm_non-ATCG_genomes.ma, sites that do not satisfy the simplex base are still excluded
        if len(uniq_nc) > 1 and tmp_all > f_base:
            cur = 0
            for key in uniq_nc.keys():
                if uniq_nc[key] > 1:
                    cur += 1

            if cur > 1 and seq[width] != '-':
                max_dic = max(uniq_nc, key=uniq_nc.get)
                pos_freq = 1 - (uniq_nc[max_dic] / tmp_all)
                tmp_freq.append(width)
                tmp_freq.append(pos_freq)
                tmp_freq.append(tmp_all)
                freqs.append(tmp_freq)
    return freqs


# rm_pis is a parameter set according to remove_fq. The parsimonious information loci above the remove_fq frequency contain only A,T,C,G sequences
def rm_N_list(alignments, rm_pis, filter_seq):
    sig_nc = {'a', 't', 'g', 'c', 'A', 'T', 'C', 'G'}
    for alignment in tqdm(alignments, desc='generating ePIS only with A,T,C,G'):
        flag = 0
        for site in rm_pis:
            if alignment.seq[site] not in sig_nc:
                flag = 1
                break
        if flag:
            continue
        else:
            SeqIO.write(alignment, filter_seq, "fasta")


# Restore the position after multiple sequence alignment
def get_ref_pos(ref_seq):
    fa_pos = []
    nt_pos = 0
    for nt in range(0, len(ref_seq)):
        if ref_seq[nt] != '-':
            # The actual sequence position is the subscript +1
            nt_pos += 1
            fa_pos.append(nt_pos)
        else:
            fa_pos.append(nt_pos)
    return fa_pos


# Output fasta files containing only parsimonious information loci
def write_pis(alignment, all_pis, pos_fa):
    pos_boolean = np.zeros(len(alignment[0]))
    pos_boolean[all_pis] = 1
    pos_boolean = list(pos_boolean)
    for n in tqdm(range(0, len(alignment)), desc='Write ePIS into fasta files'):
        pos_seq = ''
        old_seq = alignment[n].seq
        for i in range(0, len(old_seq)):
            if pos_boolean[i]:
                pos_seq = pos_seq + old_seq[i]

        pos_fa.write(">"+alignment[n].id)
        pos_fa.write("\n")
        pos_fa.write(pos_seq)
        pos_fa.write("\n")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inpath', dest='inpath', help="input file path", required=True)
    parser.add_argument('-m', '--ma', dest='ma', help="input multiple alignment result file", required=True)
    parser.add_argument('-b', '--base', dest='base', help="Number of effective bases, none if unknown", required=True)
    parser.add_argument('-r', '--remove', dest='remove', help="pis frequence cutoff of remove genomes, none if unknown, 0 if remove all non-A,T,C,G bases", required=True)
    parser.add_argument('-f', '--reference', dest='reference', help="id of reference sequence in the multiple alignment result file", required=True)
    args = parser.parse_args()

    if args.inpath[-1] != '/':
        inpath = args.inpath + "/"
    else:
        inpath = args.inpath
    try:
        align = AlignIO.read(inpath + args.ma, "fasta")
    except ValueError:
        print("Not a Multiple Sequence Alignment result file. ")
        print("Please check the input file format")
        sys.exit()

    freq_out = open(inpath + 'freq_all.txt', 'w')
    pos_out = open(inpath + 'pi_pos_all.fasta', 'w')
)
    if args.base == 'none':
        base = int(len(align)*0.8)
    else:
        base = int(args.base)
    if args.remove == 'none':
        remove = 1/pow(10, len(str(len(align)))-2)
    else:
        remove = float(args.remove)

    ref_seq = ''
    ref_flag = 0
    ref_id = args.reference
    for gn in align:
        if gn.id == ref_id:
            ref_seq = gn.seq
            ref_flag = 1
            break
    if not ref_flag:
        sys.exit("No reference genome, please check the input file")
    ref_pos = get_ref_pos(ref_seq)

    all_sites = []
    rm_sites = []
    if os.path.exists(inpath + "rm_non-ATCG_genomes.ma") and os.path.getsize(inpath + "rm_non-ATCG_genomes.ma"):
        alignments_filtered = AlignIO.read(inpath + "rm_non-ATCG_genomes.ma", "fasta")
        all_sites, rm_sites = find_ps_sites(align, ref_seq, remove, base)
        alignments_filtered = AlignIO.read(inpath + "rm_non-ATCG_genomes.ma", "fasta")

    else:
        all_sites, rm_sites = find_ps_sites(align, ref_seq, remove, base)
        filter_file = open(inpath + "rm_non-ATCG_genomes.ma", "w")
        rm_N_list(align, rm_sites, filter_file)
        filter_file.close()
        alignments_filtered = AlignIO.read(inpath + "rm_non-ATCG_genomes.ma", "fasta")

    base_filtered = int(len(alignments_filtered)*0.8)
    freq_out.write("reference")
    freq_out.write("\t")
    freq_out.write("alignment_pos")
    freq_out.write("\t")
    freq_out.write("frequence")
    freq_out.write("\t")
    freq_out.write("depth")
    freq_out.write("\n")

    # Directly generated sites start at 0. The sites subscript from the freq_all file starts at 1
    write_sites = []
    all_freq = pis_freq(alignments_filtered, all_sites, ref_seq, base_filtered)
    for freq_cur in range(0, len(all_freq)):
        tmp_index = all_freq[freq_cur]
        write_sites.append(tmp_index[0])
        orgin_pos = ref_pos[tmp_index[0]]
        freq_out.write("%d" % orgin_pos)
        freq_out.write("\t")
        align_pos = tmp_index[0] + 1
        freq_out.write("%d" % align_pos)
        freq_out.write("\t")
        freq_out.write("%.4f" % tmp_index[1])
        freq_out.write("\t")
        freq_out.write("%d" % tmp_index[2])
        freq_out.write("\n")
    write_pis(alignments_filtered, write_sites, pos_out)
    freq_out.close()
    pos_out.close()
