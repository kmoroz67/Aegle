import pandas as pd
import numpy as np
import os
from pathlib import Path
import subprocess
from itertools import islice
import shutil
import logging
import time
import argparse
import getpass
import sys
import functions
import geo_downloader
import tcga_downloader

if __name__ == "__main__":
    
    #create argparse logger
    functions.SetLogger(logger_name='argparse')

    #argparse main parameters
    logger = logging.getLogger('argparse')
    user_name = getpass.getuser()
    n_cores_upper = int(subprocess.check_output("nproc", shell = True).split()[0])-1
    available_users = subprocess.check_output("cat /etc/passwd | grep /home | cut -d: -f1", shell = True)
    available_users = available_users.decode("utf-8").strip().split()
    
    parser = argparse.ArgumentParser(description = 'TCR/BCR calculations', epilog = 'Enjoy the program! :)',
                                    formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-u', '--user', type = str, 
                help = 'Name of the user for saving intermediate results. Default: {user}.'.format(user = user_name),
                default = '{user}'.format(user = user_name))
    parser.add_argument('--dbase', type = str,
                       help = 'Name of database where data will be downloaded. Available databases: geo, tcga.', required=True)
    parser.add_argument('--n_start', type = int, help = 'Start of the range of CPUs for MiXCR calculations.', required=True)
    parser.add_argument('--n_end', type = int, help = 'End of the range of CPUs for MiXCR calculations.', required=True)
    args = parser.parse_args()
    if args.user not in available_users:
        logger.info("User {user} doesn't exist!".format(user = args.user))
        sys.exit()
    if args.dbase not in ['geo', 'tcga']:
        logger.info("Oops...Unknown database :(")
        logger.info("Available databases: geo, tcga")
        sys.exit()
    if (args.n_start > n_cores_upper or args.n_start < 0) or (args.n_end > n_cores_upper or args.n_end < 0): 
        logger.info("Invalid range for CPUs!")
        logger.info("Try to choose values between 0 and {n_cores_upper}.".format(n_cores_upper = n_cores_upper))
        sys.exit()
    if args.n_start >= args.n_end:
        logger.info("[n_start] should be less than [n_end] !")
        logger.info("Use command: python3.7 {program_name} [-h] to get more information.".format(program_name = parser.prog))
        sys.exit()
        
    #main part
    if args.dbase == 'geo':
        geo_downloader.CalculateGEO(user=args.user, 
                                    n_start=args.n_start,
                                    n_end=args.n_end)
    elif args.dbase == 'tcga':
        tcga_downloader.CalculateTCGA(user=args.user,
                                      n_start=args.n_start,
                                      n_end=args.n_end)
        
    