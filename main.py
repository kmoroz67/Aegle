import argparse
import getpass
import logging
import subprocess
import sys
from argparse import Namespace

from functions import SetLogger
import geo_downloader
import tcga_downloader


def parse_parameters() -> Namespace:

    # create argparse logger
    SetLogger(logger_name='argparse')

    # Getting the user name, the number of cores, and the available users.
    logger = logging.getLogger('argparse')
    user_name = getpass.getuser()
    n_cores_upper = int(subprocess.check_output("nproc", shell=True).split()[0]) - 1
    available_users = subprocess.check_output("cat /etc/passwd | grep /home | cut -d: -f1", shell=True)
    available_users = available_users.decode("utf-8").strip().split()

    parser = argparse.ArgumentParser(description='TCR/BCR calculations', epilog='Enjoy the program! :)',
                                     formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-u', '--user', type=str,
                        help=f'Name of the user for saving intermediate results. Default: {user_name}.',
                        default=f'{user_name}')
    parser.add_argument('--dbase', type=str,
                        help='Name of database where data will be downloaded. Available databases: geo, tcga.',
                        required=True)
    parser.add_argument('--n_start', type=int, help='Start of the range of CPUs for MiXCR calculations.', required=True)
    parser.add_argument('--n_end', type=int, help='End of the range of CPUs for MiXCR calculations.', required=True)
    args = parser.parse_args()

    if args.user not in available_users:
        logger.info(f"User {args.user} doesn't exist!")
        sys.exit()
    if args.dbase not in ['geo', 'tcga']:
        logger.info("Oops...Unknown database :(")
        logger.info("Available databases: geo, tcga")
        sys.exit()
    if (args.n_start > n_cores_upper or args.n_start < 0) or (args.n_end > n_cores_upper or args.n_end < 0):
        logger.info("Invalid range for CPUs!")
        logger.info(f"Try to choose values between 0 and {n_cores_upper}.")
        sys.exit()
    if args.n_start >= args.n_end:
        logger.info("[n_start] should be less than [n_end] !")
        logger.info(
            f"Use command: python3.7 {parser.prog} [-h] to get more information.")
        sys.exit()

    return args


def main() -> None:
    """
    It takes in a database name and a user name, and then downloads the data from the database
    """
    args = parse_parameters()
    # main part
    if args.dbase == 'geo':
        geo_downloader.CalculateGEO(user=args.user,
                                    n_start=args.n_start,
                                    n_end=args.n_end)
    elif args.dbase == 'tcga':
        tcga_downloader.CalculateTCGA(user=args.user,
                                      n_start=args.n_start,
                                      n_end=args.n_end)


if __name__ == "__main__":
    main()
