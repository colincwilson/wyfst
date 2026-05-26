# Report symbols in each column ('field') of a data file.
import sys
import configargparse
import polars as pl
from collections import Counter


def get_symbols(data_file, delim, field=None, sep=None, has_header=True):
    df = pl.read_csv(data_file, separator=delim, has_header=has_header)
    fields = list(df.columns)
    if field is not None:
        fields = [field]
    for field in fields:
        syms = Counter()
        if field != 'Morph':
            for x in df[field].to_list():
                if sep and sep != '':
                    x = x.split(sep)
                syms.update(x)
        else:
            for x in df[field].to_list():
                syms.update(x.split(';'))
        print(f'\nSymbols in field {field}')
        print(syms)
    print()


# # # # # # # # # #

if __name__ == "__main__":
    parser = configargparse.ArgumentParser(
        formatter_class=configargparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument( \
        '--file',
        required=True,
        type=str,
        help='Full path to data file.')
    parser.add_argument( \
        '--delim',
        type=str,
        default=',',
        help='Column delimiter in data file.')
    parser.add_argument( \
        '--field',
        type=str,
        help='Name of column to process (leave unspecified for all columns).')
    parser.add_argument( \
        '--sep',
        type=str,
        default='',
        help='Symbol separator within column (default empty).')
    parser.add_argument( \
        '--no_header',
        action='store_true',
        help='Data file does not have a header row.')

    args = parser.parse_args()
    has_header = not args.no_header

    get_symbols(args.file, '\t', args.field, args.sep, has_header)
