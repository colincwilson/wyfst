# Get symbols and their frequences that appear in one
# or each column (field) of a data file.
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
        syms = dict(syms)
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
        help='Symbol separator within column.')
    parser.add_argument( \
        '--header',
        action='store_true',
        help='Data file has a header row.')

    args = parser.parse_args()
    if args.delim == '\\t':
        args.delim = '\t'

    get_symbols(args.file, args.delim, args.field, args.sep, args.header)
