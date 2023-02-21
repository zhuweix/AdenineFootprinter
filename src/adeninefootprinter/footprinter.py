import argparse
from src.adeninefootprinter import *

parser = argparse.ArgumentParser(
    prog='AdenineFootprinter',
    description=("Identify Nucleosome Footprints and Accessible regions "
                 "using methylated adenine"))

subparser = parser.add_subparsers(dest='command')

indexref = subparser.add_parser('index')
predict = subparser.add_parser('predict', help='Predict Nucleosome Footprint and Accessible Region')
bamtobed = subparser.add_parser('bamtobed',
                                help='Extract Single-Molecule Footprint and Accessible Region from BAM file')
bedtobedgraph = subparser.add_parser('bedtobedgraph',
                                     help='Transform exported bed files to bedgraph')
train = subparser.add_parser('train', help='train model using methylated gDNA control')

indexref.add_argument('-r', type=str, help='Reference File in FASTA Format', required=True)
indexref.add_argument('-o', type=str, help='Name of the index file', required=True)

predict.add_argument('--bed', type=str, help='BED file with single-molecule m6A sites', required=True)
predict.add_argument('--output', type=str, help='Output bam file', required=True)
predict.add_argument('--window', type=int, help='[Optional] Window size. Default=25', default=25)
predict.add_argument('--model', type=str, help='[Optional] Model File',
                     default='../asset/default_model.data')

bamtobed.add_argument('--bam')