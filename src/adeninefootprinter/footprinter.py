import argparse

from .indexreference import index_reference
from .trainmodel import train_model
from .predictfootprint import predict_footprint
from .igv_bam_to_bed import igv_bam_to_bed
from .bed_to_bedgraph import bed_to_bedgraph


def main():
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

    predict.add_argument('--bam', type=str, help='[Recommended] BAM file with single-molecule m6A sites',
                         default='')
    predict.add_argument('--bed', type=str, help='[File1,File2,...]BED file with single-molecule m6A sites',
                         default='')
    predict.add_argument('-p', '--prefix', type=str, help='Prefix of Output BAM file', required=True)
    predict.add_argument('-m', '--model', type=str, help='[Optional] Model File',
                         default='../asset/SacCer3_EcoGII.model.json')
    predict.add_argument('--ref', type=str, help='Index File generated by footprinter index', required=True)
    predict.add_argument('--max_amb_fuse', type=int, help=('[Optional] Maximal Length to Fuse Ambiguous Region to '
                                                           'its adjacent nucleosome, Default=150. '
                                                           'Use 0 for no fusion'), default=150)
    predict.add_argument('--header', type=str, help=('[Optional] Header Info for output BAM File (ONLY for BED input)'
                                                     'BAM File Header Dict (see pysam documentation) stored by pickle. '
                                                     'Automated header will be generated if this file is not provided.'
                                                     ), default=None)
    predict.add_argument('--filter_acc', type=str, help=('[Optional] Only Export Reads with at least one '
                                                         'high-confidence Accessible Region. Default=True'
                                                         ), default=True)

    bamtobed.add_argument('-b', '--bam', type=str, help='BAM file from prediction', required=True)
    bamtobed.add_argument('-f', '--footprint', type=str, help='Output BED file for nucleosome footprint', required=True)
    bamtobed.add_argument('-a', '--accessible', type=str, help='Output BED file for accessible region', required=True)

    bedtobedgraph.add_argument('-b', '--bed', help='[File1,File2,...] BED file(s) transformed from BAM prediction file', required=True)
    bedtobedgraph.add_argument('-o', '--output', type=str, help='Output BEDGRAPH file for density', required=True)
    bedtobedgraph.add_argument('-r', '--ref', type=str, help='Index File generated by footprinter index', required=True)
    bedtobedgraph.add_argument('-n', '--name', type=str, help='[Optional] Name of the Track. Default=Feature',
                               default='Feature')
    train.add_argument('-b', '--bed', type=str, help='[File1,File2,...] BED file(s) for m6A sites in gDNA control',
                       required=True)
    train.add_argument('-w', '--window', type=int, help='Window size. Default=25, Minimal=5', default=25)
    train.add_argument('-p', '--prefix', type=str, help=('Prefix for Output Files. '
                                                         'Model File in [prefix].json. '
                                                         'Raw local m6A info in [prefix].Window=[window].csv'),
                       required=True)
    train.add_argument('--fp', type=float, help='False Positive Rate for prediction. Default=0.01', default=0.01)
    train.add_argument('--ref', type=str, help='Index File generated by footprinter index', required=True)

    args = parser.parse_args()
    if args.command == 'index':
        ref = args.r
        output = args.o
        print('Making Reference Dictionary for {}'.format(ref))
        index_reference(ref=ref, output=output)
        print('Finished')
    elif args.command == 'train':
        bed = args.bed
        window = args.window
        prefix = args.prefix
        fprate = args.fp
        ref = args.ref
        train_model(
            bed_list=bed,
            window=window,
            prefix=prefix,
            fprate=fprate,
            ref=ref
        )
    elif args.command == 'predict':
        bam = args.bam
        bed = args.bed
        ref = args.ref
        model = args.model
        max_amb_fuse = args.max_amb_fuse
        header = args.header
        filter_acc = args.filter_acc
        prefix = args.prefix
        predict_footprint(
            bam=bam,
            bed=bed,
            ref=ref,
            model=model,
            prefix=prefix,
            max_amb_fuse=max_amb_fuse,
            header_file=header,
            filter_acc=filter_acc)
    elif args.command == 'bamtobed':
        bam = args.bam
        footprint = args.footprint
        accessible = args.accessible
        igv_bam_to_bed(
            bam=bam,
            acc=accessible,
            nuc=footprint
        )
    elif args.command == 'bedtobedgraph':
        bed = args.bed
        output = args.output
        ref = args.ref
        name = args.name
        bed_to_bedgraph(
            bed=bed,
            bedgraph=output,
            name=name,
            ref=ref
        )


if __name__ == '__main__':
    main()