#!/usr/bin/env python
# coding=utf-8
# __author__ = 'Yunchao Ling'

from Bio import SeqIO
from tqdm import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import click
import shutil
import os
import multiprocessing
import sys
import re


@click.command()
@click.argument("fasta_file", type=click.Path(exists=True))
@click.argument("ref_id", nargs=1)
@click.argument("out_dir", nargs=1)
def multi_mafft(fasta_file: str, ref_id: str, out_dir: str):
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.mkdir(out_dir)

    # print("拆分出MAFFT输入文件")
    print("Split out the MAFFT input file ")
    mafft_input_dir = os.path.join(out_dir, "mafft_input")
    os.mkdir(mafft_input_dir)
    gen_mafft_input(fasta_file, ref_id, mafft_input_dir)

    # print("执行多进程并行MAFFT")
    print("Perform multi-process parallel MAFFT ")
    mafft_output_dir = os.path.join(out_dir, "mafft_output")
    os.mkdir(mafft_output_dir)
    p = multiprocessing.Pool(multiprocessing.cpu_count())
    for fl in tqdm(os.listdir(mafft_input_dir), desc="Virus"):
        fl_path = os.path.join(mafft_input_dir, fl)
        if os.path.isfile(fl_path) and fl.endswith(".fasta"):
            out_name = fl[:fl.rindex(".")] + ".ma"
            p.apply_async(do_mafft, args=(fl_path, os.path.join(mafft_output_dir, out_name)))
    p.close()
    p.join()

    # print("合并为比对完成的MA文件")
    print("Merge into a matched MA file ")
    ma_file_name = os.path.basename(fasta_file)
    ma_file_name = ma_file_name[:ma_file_name.rindex(".")] + ".ma"
    seqrecords = []
    for ma in tqdm(os.listdir(mafft_output_dir), desc="Virus"):
        if ma.endswith(".ma"):
            tqdm.write(ma)
            seq_id, seq = get_ma(os.path.join(mafft_output_dir, ma))
            seqrecord = SeqRecord(Seq(seq, IUPAC.IUPACAmbiguousDNA), id=seq_id, name="", description="")
            seqrecords.append(seqrecord)
    SeqIO.write(seqrecords, os.path.join(out_dir, ma_file_name), "fasta")

    # print("删除中间临时文件")
    print("Delete intermediate temporary files")
    shutil.rmtree(mafft_input_dir)
    shutil.rmtree(mafft_output_dir)


def gen_mafft_input(fasta_file: str, ref_id: str, mafft_input_dir: str):
    ref_seqrecord = None
    seqrecords = SeqIO.parse(fasta_file, "fasta")
    for seqrecord in seqrecords:
        if seqrecord.id == ref_id:
            ref_seqrecord = seqrecord
            break

    if ref_seqrecord is None:
        # print("Fasta文件中不包含reference: %s, 程序退出" % ref_id)
        print("No reference in Fasta file: %s, Program exit" % ref_id)
        sys.exit(1)
    else:
        for seqrecord in tqdm(SeqIO.parse(fasta_file, "fasta"), desc="Virus"):
            seq_id = seqrecord.id
            seq = str(seqrecord.seq)
            seq = re.sub(r'\.', "-", seq)
            seqrecord.seq = Seq(seq, IUPAC.IUPACAmbiguousDNA)
            SeqIO.write([ref_seqrecord, seqrecord], os.path.join(mafft_input_dir, seq_id + ".fasta"), "fasta")


def do_mafft(infile: str, out_file: str):
    os.system("mafft --thread 1 --quiet " + infile + " > " + out_file)


def get_ma(ma_file):
    seqrecords = list(SeqIO.parse(ma_file, "fasta"))
    seq_id = seqrecords[1].id
    ref_seq = str(seqrecords[0].seq)
    fa_seq = str(seqrecords[1].seq)
    ma_seq = ""
    for i in range(len(ref_seq)):
        if ref_seq[i] != "-":
            ma_seq += fa_seq[i]
    return seq_id, ma_seq


if __name__ == '__main__':
    multi_mafft()
