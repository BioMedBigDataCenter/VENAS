# __title__   = Find parsimony informative sites
# __created__ = 2020-09-09
# __author__  = Ruifang Cao
# __python__  = 3.+

"""
查找alignment结果中的简并信息位点位点所在位置，
这个位点至少有两种以上碱基的情况，同时不同碱基的出现评率大于2
v1.6:将remove非“ATGC”的碱基的基因组去除脚本和简约信息位点脚本合并为一个脚本，只查找一遍PIS，
先根据PIS筛选掉一些基因组，然后计算这些PIS对应的评率
这一版本的文件会比pis_remove_N.py剔除的基因组更多，因为数据pis_remove_N.py脚本仅剔除对应频率在是N或者-的碱基
而这个脚本简并碱基，例如Y，S的基因组也剔除
"""

from Bio import AlignIO, SeqIO
import os
from tqdm import tqdm
import numpy as np
import argparse
import sys


# 第一轮查找，all_pis为所有有效简约信息文件，rm_pis是根据用户给的阈值筛选一定频率以上用作晒酸基因组在该pis位点均是ATGC碱基的pis位点
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
        # 有效碱基数大于base
        if len(uniq_nc) > 1 and tmp_all > epis_base:
            cur = 0
            for key in uniq_nc.keys():
                if uniq_nc[key] > 1:
                    cur += 1

            if cur > 1 and seq[width] != '-':
                max_dic = max(uniq_nc, key=uniq_nc.get)
                pos_freq = 1 - (uniq_nc[max_dic] / tmp_all)
                all_pis.append(width)
                # remove_fq 是用户定义的需要删除那些频率以内的碱基是非 “ATCG”的情况
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
        # fasta文件前n条序列是outgroup序列，从outgroup序列后统计简约位点
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
        # 有效碱基数大于1000 and tmp_all > 1000
        # rm_non-ATCG_genomes.ma中不满足简并碱基的点仍排除
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


# rm_pis根据remove_fq设置了参数，remove_fq频率以上的简约信息位点只包含A,T,C,G序列
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


# 还原多序列比对后的位置
def get_ref_pos(ref_seq):
    fa_pos = []
    nt_pos = 0
    for nt in range(0, len(ref_seq)):
        if ref_seq[nt] != '-':
            # 实际位置是下标加1
            nt_pos += 1
            fa_pos.append(nt_pos)
        else:
            fa_pos.append(nt_pos)
    return fa_pos


# 根据简并信息位点输出新的fasta文件，出简并信息位点位点外的其他位点均为-
def write_pis(alignment, all_pis, pos_fa):
    #生成一个长度和alignment结果各条序列长度一致，且sites位点为1。其他位点为0的向量
    pos_boolean = np.zeros(len(alignment[0]))
    pos_boolean[all_pis] = 1
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
    # 参数设置
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inpath', dest='inpath', help="input file path", required=True)
    parser.add_argument('-m', '--ma', dest='ma', help="input multiple alignment result file", required=True)
    parser.add_argument('-b', '--base', dest='base', help="Number of effective bases, none if unknown", required=True)
    parser.add_argument('-r', '--remove', dest='remove', help="pis frequence cutoff of remove genomes, none if unknown, 0 if remove all non-A,T,C,G bases", required=True)
    parser.add_argument('-f', '--reference', dest='reference', help="id of reference sequence in the multiple alignment result file", required=True)
    args = parser.parse_args()
    # 加载参数
    # 设置输入路径
    if args.inpath[-1] != '/':
        inpath = args.inpath + "/"
    else:
        inpath = args.inpath
    # 读取输入aligment文件
    try:
        align = AlignIO.read(inpath + args.ma, "fasta")
    except ValueError:
        print("Not a Multiple Sequence Alignment result file. ")
        print("Please check the input file format")
        sys.exit()

    # 设置数据freq和pi_pos文件名
    freq_out = open(inpath + 'freq_all.txt', 'w')
    pos_out = open(inpath + 'pi_pos_all.fasta', 'w')

    # 设置有效碱基数 (判断用户输入的必须是整数)
    if args.base == 'none':
        base = int(len(align)*0.8)
    else:
        base = int(args.base)
    # 设置需要排除ATCG情况的pis的频率阈值 (判断用户输入的必须是float)
    if args.remove == 'none':
        remove = 1/pow(10, len(str(len(align)))-2)
    else:
        remove = float(args.remove)

    ref_seq = ''
    ref_flag = 0
    ref_id = args.reference
    for gn in align:
        # OEAV139851, AY390556.1, OEAV118620，OEAV121000
        if gn.id == ref_id:
            ref_seq = gn.seq
            ref_flag = 1
            break
    if not ref_flag:
        sys.exit("No reference genome, please check the input file")
    ref_pos = get_ref_pos(ref_seq)

    # 根据阈值去除pis中不是ATCG的序列
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

    # 加载去除非A,T,C,G后的序列后，根据过滤后的alignment文件包含基因组的个数设置对应的base值。
    base_filtered = int(len(alignments_filtered)*0.8)

    # 输出频率文件，第一列为在reference中的位置，第二列为alignment结果位置，第三列为位点频率
    freq_out.write("reference")
    freq_out.write("\t")
    freq_out.write("alignment_pos")
    freq_out.write("\t")
    freq_out.write("frequence")
    freq_out.write("\t")
    freq_out.write("depth")
    freq_out.write("\n")

    # 这里需要注意直接生成的sites是从0开始，从freq_all文件得到的sites下标是从1开始
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
