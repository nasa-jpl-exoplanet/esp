'''monitors the pipeline and corrects the system when it sees problems'''

import argparse
import os
import sys

from urllib.parse import urlparse

from . import fsm
from . import mia


def _file(arg: str):
    if arg is not None:
        if not os.path.isfile(arg):
            raise argparse.ArgumentTypeError(
                f'Path "{arg}" is not a valid file'
            )
    return arg


def _url(arg: str):
    try:
        result = urlparse(arg)
        if all((result.scheme, result.netloc)):
            return arg
        raise argparse.ArgumentTypeError(
            f"Invalid URL: '{arg}' does not have a scheme or network location."
        )
    except ValueError as ve:
        raise argparse.ArgumentTypeError(f"Invalid URL format: '{arg}'") from ve


def tyrannize():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_arugment(
        '-ca',
        '--ca-path',
        default=None,
        required=False,
        type=_file,
        help='path to the CA cerficates file for this machine',
    )
    ap.add_argument(
        '-e',
        '--emails',
        required=True,
        type=str,
        help='email addresses to notify in a , delimited string',
    )
    ap.add_argument(
        '-u',
        '--url',
        default='https://excalibur.jpl.nasa.gov:8080',
        type=_url,
        help='the URL for the pipeline [%(default)s]',
    )
    sub = ap.add_subparsers(dest="command", help='Available Commands')
    sap = sub.add_parser('fsm', help=fsm.__doc__)
    sap.add_argument(
        '-t',
        '--threshold',
        required=True,
        type=int,
        help='the duration in seconds that the FSM is allowed to be stuck in a state other than "running" before restarting the pipeline',
    )
    sap.set_defaults(func=fsm.check)
    nap = sub.add_parser('mia-farm', help=mia.farm.__doc__)
    nap.add_argument(
        '-n',
        '--nodes',
        nargs='+',
        required=True,
        type=int,
        help='mentor node numbers',
    )
    nap.set_defaults(func=mia.farm)
    wap = sub.add_parser('mia-worker', help=mia.worker.__doc__)
    wap.add_argument(
        '-n',
        '--nodes',
        nargs='+',
        required=True,
        type=int,
        help='mentor node numbers',
    )
    wap.add_argument(
        '-r',
        '--replicas',
        required=True,
        type=int,
        help='the number of workers on each node',
    )
    wap.set_defaults(func=mia.worker)
    args = ap.parse_args()
    if hasattr(args, 'func'):
        sys.exit(args.func(args))
    else:
        ap.print_help()
        sys.exit(-1)


if __name__ == '__main__':
    tyrannize()
