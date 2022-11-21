import sys
import argparse
from . import _program, crop_and_filter_plate, pixel_counts, make_dir
from clint.textui import puts_err, indent, colored
import os
import sys


def main(args = sys.argv[1:]):
    parser = argparse.ArgumentParser(prog = _program)

    parser.add_argument('images', type=argparse.FileType('r'), nargs='+')

    parser.add_argument("--radius",
                        help="Set the radius (px) of the plate",
                        type=int,
                        default=930)

    parser.add_argument("-e",
                        "--crop",
                        help="Additional crop from edge",
                        type=int,
                        default=100)

    parser.add_argument("-s",
                        "--small",
                        help="Filter out particles smaller than specified",
                        type=int,
                        default=100)

    parser.add_argument("-l",
                        "--large",
                        help="Filter out particles larger than specified",
                        type=int,
                        default=1500)

    parser.add_argument("-r",
                    "--rev",
                    help="Reverse the sign of the CI index",
                    action="store_true")

    parser.add_argument("--header",
                        help="Output a header",
                        default=False,
                        action="store_true")

    parser.add_argument("-c",
                    "--center",
                    help="Divisor for defining center",
                    type=int,
                    default=5)

    parser.add_argument("-d",
                        "--debug",
                        help="Output debug information",
                        action="store_true",
                        default=False)

    args = parser.parse_args(args)

    if args.header:
        print("fname\tq1\tq2\tq3\tq4\tn\ttotal_q\ttotal\tci")

    radius_range = [args.radius - 20, args.radius + 20]

    if args.debug:
        make_dir("debug")

    for img in args.images:
        fname = os.path.splitext(os.path.basename(img.name))[0]
        with indent(4):
            puts_err(colored.blue("\nProcessing " + img.name))
            img = crop_and_filter_plate(img.name,
                                        radius_range,
                                        small=args.small,
                                        large=args.large,
                                        extra_crop=args.crop,
                                        debug=args.debug)

            puts_err(colored.blue("Use pixel method for calculation"))
            result = pixel_counts(img, args.center)
            if args.rev:
                result[-1] = result[-1]*-1.0
            print(fname + "\t" + '\t'.join(map(str, result)))
            sys.stdout.flush()



if __name__ == '__main__':
    main()
