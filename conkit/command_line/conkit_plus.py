import argparse
import conkit.plot
import conkit.io

def parse_arguments():
    parser = argparse.ArgumentParser(description='ConKit: Interactive contact plot with prediction data visualisation',
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("confile", help="Path to contact file", type=str)
    parser.add_argument("conformat", help="Format of contact file", type=str)

    parser.add_argument("seqfile", help="Path to sequence file", type=str)
    parser.add_argument("seqformat", help="Format of sequence file", type=str)

    parser.add_argument("-reffile", help="Path to pdb reference file", nargs="?", default=None, type=str)
    parser.add_argument("-refformat", help="Format of reference file", nargs="?", default=None, type=str)

    parser.add_argument("-ssfile", help="Path to ss file", nargs="?", default=None, type=str)
    parser.add_argument("-ssformat", help="Format of ss file", nargs="?", default=None, type=str)

    parser.add_argument("-memfile", help="Path to mem file", nargs="?", default=None, type=str)
    parser.add_argument("-memformat", help="Format of mem file", nargs="?", default=None, type=str)

    parser.add_argument("-conservationfile", help="Path to conservation file", nargs="?", default=None, type=str)
    parser.add_argument("-conservationformat", help="Format of conservation file", nargs="?", default=None, type=str)

    parser.add_argument("-disorderfile", help="Path to disorder file", nargs="?", default=None, type=str)
    parser.add_argument("-disorderformat", help="Format of disorder file", nargs="?", default=None, type=str)

    args = parser.parse_args()

    return args


def main():
    args = parse_arguments()

    seq = conkit.io.read(args.seqfile, args.seqformat).top
    conpred = conkit.io.read(args.confile, args.conformat).top
    test_plot = conkit.plot.ContactPlotPlusFigure(seq=seq, conpred=conpred)

    if args.refformat is not None:
        reference = conkit.io.read(args.reffile, args.refformat).top
        test_plot.add_reference_layer(reference)

    if args.ssformat is not None:
        sspred = conkit.io.read(args.ssfile, args.ssformat)
        test_plot.add_sspred_layer(sspred)

    if args.conservationformat is not None:
        conservpred = conkit.io.read(args.conservationfile, args.conservationformat)
        test_plot.add_conservpred_layer(conservpred)

    if args.memformat is not None:
        mempred = conkit.io.read(args.memfile, args.memformat)
        test_plot.add_mempred_layer(mempred)

    if args.disorderformat is not None:
        dispred = conkit.io.read(args.disorderfile, args.disorderformat)
        test_plot.add_dispred_layer(dispred)

    test_plot.plot(fname='/Users/shahrammesdaghi/Downloads/test.html')

if __name__ == "__main__":

    main()

