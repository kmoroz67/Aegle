import pandas as pd
import numpy as np
import os
from pathlib import Path
import itertools
import subprocess
from itertools import islice
import shutil
import logging
import time
import argparse
import getpass
import sys


def SetLogger(logger_name):
    """
    Create custom logger and set its configuration
    
    :param logger_name: name of created logger
    """
     
    #Create a custom logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    
    #Create handlers
    c_handler = logging.StreamHandler()
    f_handler = logging.FileHandler('/uftp/Blood/db_calc_pipeline/tmp/pipeline.log')
    c_handler.setLevel(logging.INFO)
    f_handler.setLevel(logging.WARNING)
    
    #Create formatters and add them to handlers
    if logger_name == 'argparse': 
        c_format = logging.Formatter('%(message)s')
    else:
        c_format = logging.Formatter('%(levelname)s: %(message)s')
    f_format = logging.Formatter('%(asctime)s::%(name)s::%(levelname)s::%(message)s',  datefmt='%Y-%m-%d %H:%M:%S')
    c_handler.setFormatter(c_format)
    f_handler.setFormatter(f_format)
    
    #Add handlers to the logger
    logger.addHandler(c_handler)
    logger.addHandler(f_handler)

def DownloadFasterq(run, sample_id, dataset, run_path, df, df_stats, user, layout, logger_name='fasterq'):
    """
    Download .fastq files from geo, unzip it and make its validation

    :param run: name of run (SRR from GEO) in calculated dataset
    :param sample_id: name of sample (SRX from GEO) in calculated dataset
    :param dataset: name of dataset 
    :param run_path: path to SRR folder (where the data will be downloaded)
    :param df: dataframe of calculation list (list of datasets to download)
    :param df_stats: cumulative dataframe with statistics for one sample
    :param user: user name to get access to fasterq-dump tool
    :param layout: type of run (paired or single)
    :param logger_name: name of created logger (by default logger_name = "fasterq")
    :return: df_stats - cumulative dataframe with statistics for one sample
    """
    
    
    logger = logging.getLogger(logger_name)
    logger.info('start fasterq-dump for {run}'.format(run=run))
    
    n_reps = df.n_repeats[df.Run == run].values[0]  # number of attempts to download this geo data
    counter = 0  # counter of download repeats 
    
    while True:
        try:
            template = '/home/{user}/ncbi/usr/local/ncbi/sra-tools/bin/fasterq-dump-orig.2.10.9 {run} -O {path_to_dir} -t {path_to_dir} -f'
            logger.info(template.format(user=user, path_to_dir=str(run_path), run=run))
            subprocess.run(template.format(user=user, path_to_dir=str(run_path), run=run), shell=True, check=True)
            break
        except:
            if counter == 25:
                mark = 'FAILED'
                df_stats = pd.concat([df_stats,
                                      pd.DataFrame(np.array([run, sample_id, dataset, mark, n_reps + 1]).reshape(1, -1),
                                                   columns=['Run', 'Sample', 'Dataset', 'Mark', 'N_repeats'])])
                subprocess.run('rm -r {}'.format(str(run_path)), shell=True)
                return df_stats
            counter += 1
            logger.warning('{} failed attempts, local run'.format(counter))
            time.sleep(60)

    mark = 'SUCCESS'
    fq_paths = {}
    if layout == 'SINGLE':
        fq_paths['SINGLE'] = [Path(run_path, run + '.fastq')]
    elif layout == 'PAIRED':
        fq_paths['PAIRED'] = [Path(run_path, run + '_{i}.fastq'.format(i=i)) for i in (1, 2)]
    if not all([fq_path.exists() for fq_path in fq_paths[layout]]):
        mark = 'FAILED'

    df_stats = pd.concat([df_stats,
                          pd.DataFrame(np.array([run, sample_id, dataset, mark, n_reps + 1]).reshape(1, -1),
                                       columns=['Run', 'Sample', 'Dataset', 'Mark', 'N_repeats'])])
    logger.info('done fasterq-dump for {run}'.format(run=run))
    return df_stats

  
def MergeFastqFiles(df_stats, sample_path, sample_id, dataset, layout, logger_name='merging'):
    """
    Merge all .fastq files (all runs) located inside one sample and make its validation. 
    Paired runs are merged into two separate files named sample_1.fastq and sample_2.fastq respectively.
    Single runs are merged into file named sample.fastq.

    :param df_stats: cumulative dataframe with statistics for one sample
    :param sample_path: path to sample folder (SRX for GEO) with downloaded data
    :param sample_id: name of sample (SRX from GEO) in calculated dataset
    :param dataset: name of dataset 
    :param layout: type of run (paired or single)
    :param logger_name: name of created logger (by default logger_name = "merging")
    """
    
    logger = logging.getLogger(logger_name)
    logger.info('start merging, local run')
    
    if layout == 'PAIRED':
        with open(sample_path / 'sample_1.fastq', 'w') as f_out_1:
            with open(sample_path / 'sample_2.fastq', 'w') as f_out_2:
                for root, dirs, files in os.walk(sample_path, topdown=True):
                    for f in [Path(root, f_name) for f_name in files if f_name.endswith('.fastq')]:
                        if str(f).endswith('_1.fastq'):
                            with open(f, 'r') as f_in_1:
                                shutil.copyfileobj(f_in_1, f_out_1)
                        elif str(f).endswith('_2.fastq'):
                            with open(f, 'r') as f_in_2:
                                shutil.copyfileobj(f_in_2, f_out_2)
        subprocess.run('chmod 777 {}'.format(str(sample_path / 'sample_1.fastq')), shell=True)
        subprocess.run('chmod 777 {}'.format(str(sample_path / 'sample_2.fastq')), shell=True)
        for item in sample_path.iterdir():
            if item.is_dir():
                subprocess.run('rm -r {}'.format(item), shell=True, check=True) 
                
    elif layout == 'SINGLE':
        with open(sample_path / 'sample.fastq', 'w') as f_out:
            for root, dirs, files in os.walk(sample_path, topdown=True):
                for f in [Path(root, f_name) for f_name in files if f_name.endswith('.fastq')]:
                    with open(f, 'r') as f_in:
                        shutil.copyfileobj(f_in, f_out)
        subprocess.run('chmod 777 {}'.format(str(sample_path / 'sample.fastq')), shell=True)
        for item in sample_path.iterdir():
            if item.is_dir():
                subprocess.run('rm -r {}'.format(item), shell=True, check=True) 
    else:
        logger.error('wrong layout {}'.format(layout))
    logger.info('merge done, local run')


def make_small_fastq(sample_path, layout, sample_id='sample', n=32000000, logger_name='optitype'):
    """
    Split merged .fastq file into separate parts and iterate through them

    :param sample_path: path to sample folder (SRX for GEO) with downloaded data
    :param layout: type of run (paired or single)
    :param sample_id: name of sample with its merged runs (by default sample_id = "sample")
    :param n: number of rows after each splitting (by default n = 16000000)
    :param logger_name: name of created logger (by default logger_name = "optitype")
    """

    logger = logging.getLogger(logger_name)
    logger.info('start fastq for optitype')
    
    output_path = Path(sample_path)
    
    if layout == 'SINGLE':
        fastq_path = output_path / '{sam}.fastq'.format(sam=sample_id)
        with fastq_path.open('r') as f:
            next_n_lines = list(islice(f, n))
            with (output_path / 'sample_optitype.fastq').open('w') as out:
                out.writelines(next_n_lines)
        subprocess.run('chmod 777 {}'.format(str(output_path / 'sample_optitype.fastq')), shell=True)
        
    elif layout == "PAIRED":
        
        fastq_path_1 = output_path / '{sam}_1.fastq'.format(sam=sample_id)
        fastq_path_2 = output_path / '{sam}_2.fastq'.format(sam=sample_id)

        with fastq_path_1.open('r') as f:
            next_n_lines = list(islice(f, n))
            with (output_path / 'sample_optitype_1.fastq').open('w') as out:
                out.writelines(next_n_lines)
        with fastq_path_2.open('r') as f:
            next_n_lines = list(islice(f, n))
            with (output_path / 'sample_optitype_2.fastq').open('w') as out:
                out.writelines(next_n_lines)
                
        subprocess.run('chmod 777 {}'.format(str(output_path / 'sample_optitype_1.fastq')), shell=True)
        subprocess.run('chmod 777 {}'.format(str(output_path / 'sample_optitype_2.fastq')), shell=True)
        
    else:
        logger.error('wrong layout {}'.format(layout))
    logger.info('done fastq for optitype')


def RunOptitype(sample_path, output_path, sample_id, dataset, layout, logger_name='optitype'):
    """
    Run Optitype
    
    :param sample_path: path to sample folder (SRX for GEO)
    :param output_path: path to file with results of optitype validation
    :param sample_id: name of sample (SRX from GEO) in calculated dataset
    :param dataset: name of dataset 
    :param layout: type of run (paired or single)
    :param logger_name: name of created logger (by default logger_name = "optitype")
    """
    
    logger = logging.getLogger(logger_name)
    logger.info('start optitype, local run')
    
    templates = {}
    templates[
        'optitype'] = 'docker run --rm --user $(id -u) -v {path_to_dir}:/data -t fred2/optitype -i sample_optitype_1.fastq sample_optitype_2.fastq -r -o /data/'

    templates['optitype-se'] = 'docker run --rm --user $(id -u) -v {path_to_dir}:/data -t fred2/optitype -i sample_optitype.fastq -r -o /data/'

    make_small_fastq(sample_path=sample_path, layout=layout, sample_id='sample', n=32000000, logger_name='optitype')
    success = True
    if layout == 'PAIRED':
        try:
            subprocess.run(templates['optitype'].format(path_to_dir=sample_path),
                           shell=True,
                           check=True)
        except:
            success = False
            logger.warning('can not calculate HLA')
            
        subprocess.run('rm {}'.format(str(Path(sample_path, 'sample_optitype_1.fastq'))), shell=True)
        subprocess.run('rm {}'.format(str(Path(sample_path, 'sample_optitype_2.fastq'))), shell=True)

    elif layout == 'SINGLE':
        try:
            subprocess.run(templates['optitype-se'].format(path_to_dir=sample_path),
                           shell=True,
                           check=True)
        except:
            success = False
            logger.warning('can not calculate HLA')
                  
        subprocess.run('rm {}'.format(str(Path(sample_path, 'sample_optitype.fastq'))), shell=True)
    else:
        success = False
        with open(output_path, 'a') as f:
            logger.error('wrong layout {}'.format(layout))
            f.write('{},{},{},failed\n'.format(sample_id, dataset, layout))
    
    # Это переименовывание папок после того как отработал optitype
    if success:  
        p = Path(sample_path)
        folders = [x for x in p.iterdir() if x.is_dir() and x.name != 'results']
        path_to_res = folders[0] / (folders[0].name + '_result.tsv')
        path_to_res.rename(folders[0] / ('hla-I.tsv'))
        path_to_pdf = folders[0] / (folders[0].name + '_coverage_plot.pdf')
        path_to_pdf.rename(folders[0] / ('coverage_plot.pdf'))
        folders[0].rename(folders[0].parent / 'hla')

    if Path(sample_path, 'hla').exists():
        with open(output_path, 'a') as f:
            f.write('{},{},{},success\n'.format(sample_id, dataset, layout))
    else:
        Path(sample_path, 'hla').mkdir(parents=True, exist_ok=True)
        with open(output_path, 'a') as f:
            logger.warning('Optitype failed for: {}'.format(sample_id))
            f.write('{},{},{},failed\n'.format(sample_id, dataset, layout))
    logger.info('optitype done, local run')


def RunMixcr(n_start, n_end, sample_path, output_path, sample_id, dataset, layout, logger_name='mixcr'):
    """
    Run MiXCR
    
    :param n_start: start of the range of CPUs for MiXCR calculations
    :param n_end: end of the range of CPUs for MiXCR calculations
    :param sample_path: path to sample folder (SRX for GEO) 
    :param output_path: path to file with results of MiXCR validation
    :param sample_id: name of sample (SRX from GEO) in calculated dataset
    :param dataset: name of dataset 
    :param layout: type of run (paired or single)
    :param logger_name: name of created logger (by default logger_name = "mixcr")
    """
        
    logger = logging.getLogger('mixcr')   
    logger.info('start mixcr')
    
    templates = {}
    templates[
        'mixcr'] = 'docker run --rm --user $(id -u) --cpuset-cpus {n_start}-{n_end} -v {path_to_dir}:/inputs {image} /inputs/{sample_id}_1.fastq /inputs/{sample_id}_2.fastq /inputs/RESULTS.tar.gz'
    templates[
        'mixcr-se'] = 'docker run --rm --user $(id -u) --cpuset-cpus {n_start}-{n_end} -v {path_to_dir}:/inputs {image} /inputs/{sample_id}.fastq hs /inputs/RESULTS.tar.gz'

    if layout == 'PAIRED':
        subprocess.run(templates['mixcr'].format(path_to_dir=sample_path,
                                                 n_start=n_start, n_end=n_end, sample_id='sample',
                                                 image='028257207274.dkr.ecr.us-east-1.amazonaws.com/hippocrates/mixcr:1.3'),
                       shell=True,
                       check=True)

    elif layout == 'SINGLE':
        subprocess.run(templates['mixcr-se'].format(path_to_dir=sample_path,
                                                    n_start=n_start, n_end=n_end, sample_id='sample',
                                                    image='028257207274.dkr.ecr.us-east-1.amazonaws.com/hippocrates/mixcr-se:1.2'),
                       shell=True,
                       check=True)
    else:
        with open(output_path, 'a') as f:
            logger.error('wrong layout {}'.format(layout))
            f.write('{},{},{},failed\n'.format(sample_id, dataset, layout))

    if Path(sample_path, 'RESULTS.tar.gz').exists():
        with open(output_path, 'a') as f:
            f.write('{},{},{},success\n'.format(sample_id, dataset, layout))
        subprocess.run('tar -zxvf {} -C {}'.format(str(Path(sample_path, 'RESULTS.tar.gz')), str(sample_path)),
                       shell=True)
        subprocess.run('mv {} {}'.format(str(Path(sample_path, 'tmp', 'results')), str(Path(sample_path))),
                       shell=True)
        subprocess.run('rm -r {}'.format(str(Path(sample_path, 'tmp'))), shell=True)
        subprocess.run('rm {}'.format(str(Path(sample_path, 'RESULTS.tar.gz'))), shell=True)
        
    else:
        with open(output_path, 'a') as f:
            logger.warning('Mixcr failed for: {}'.format(sample_id))
            f.write('{},{},{},failed\n'.format(sample_id, dataset, layout))
            
    logger.info('mixcr done')

if __name__ == "__main__":
    
    #create loggers
    loggers_used = ['argparse', 'main', 'fasterq', 
                    'merging', 'optitype', 'mixcr']
    for logger in loggers_used:
        SetLogger(logger)
    
    #argparse main parameters
    logger = logging.getLogger('argparse')
    user_name = getpass.getuser()
    n_cores_upper = int(subprocess.check_output("nproc", shell = True).split()[0])-1
    available_users = subprocess.check_output("cat /etc/passwd | grep /home | cut -d: -f1", shell = True)
    available_users = available_users.decode("utf-8").strip().split()
    
    parser = argparse.ArgumentParser(description = 'TCR/BCR calculations', epilog = 'Enjoy the program! :)')
    parser.add_argument('-u', '--user', type = str, 
                        help = 'Name of the user for saving intermediate results. Default: {user}'.format(user = user_name),
                        default = '{user}'.format(user = user_name))
    parser.add_argument('--n_start', type = int, help = 'Start of the range of CPUs for MiXCR calculations.', 
                        required = True)
    parser.add_argument('--n_end', type = int, help = 'End of the range of CPUs for MiXCR calculations.', 
                        required = True)
    args = parser.parse_args()
    if args.user not in available_users:
        logger.info("User {user} doesn't exist!".format(user = args.user))
        sys.exit()
    if (args.n_start > n_cores_upper or args.n_start < 0) or (args.n_end > n_cores_upper or args.n_end < 0): 
        logger.info("Invalid range for CPUs!")
        logger.info("Try to choose values between 0 and {n_cores_upper}.".format(n_cores_upper = n_cores_upper))
        sys.exit()
    if args.n_start >= args.n_end:
        logger.info("[n_start] should be less than [n_end] !")
        logger.info("Use command: python3.7 {program_name} [-h] to get more information.".format(program_name = parser.prog))
        sys.exit()
        
    
    #start the main cycle of list check
    logger = logging.getLogger('main')
    parent_path = Path('/uftp/Blood/db_calc_pipeline/')
    local_path = Path('/home/{user}/mixcr_calc/'.format(user=args.user))
    while True:
        df = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/ann_calc_slice.csv')
        last_SRX = df.Sample.tail(1).item()
        datasets = df.Dataset.unique()
        np.random.shuffle(datasets)
        for dataset in datasets:
            dataset_path = parent_path / dataset
            dataset_path.mkdir(parents=True, exist_ok=True)
            for sample_id in df.Sample[df.Dataset == dataset]:

                sample_path_uftp = parent_path / dataset / sample_id
                sample_path_local = local_path / dataset / sample_id

                if Path(sample_path_uftp, 'hla').exists() and Path(sample_path_uftp, 'results').exists():
                    need_to_calculate = False
                    continue
                else:  
                    df_process = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/download_process_samples.csv')
                    if sample_id in df_process['Sample'].values:
                        need_to_calculate = False
                        continue      
                    else:      
                        with open('/uftp/Blood/db_calc_pipeline/tmp/download_process_samples.csv', 'a') as f:
                            f.write('{}\n'.format(sample_id))
                        need_to_calculate = True
                        break  #####################################
                        
            if not need_to_calculate:
                continue
                                
            sample_path_local.mkdir(parents=True, exist_ok=True)
            subprocess.run('chmod 777 {}'.format(str(sample_path_local)), shell=True)

            statistics_df = pd.DataFrame(columns=['Run', 'Sample', 'Dataset', 'Mark', 'N_repeats'])
            for run in df.Run[df.Sample == sample_id]:
                layout = df.Layout[df.Run == run].values[0]
                layout = layout.upper()  # if layout will be in lowercase
                run_path = sample_path_local / run
                run_path.mkdir(parents=True, exist_ok=True)
                statistics_df = DownloadFasterq(run=run, sample_id=sample_id,
                                                dataset=dataset,
                                                run_path=run_path, df=df, 
                                                df_stats=statistics_df, 
                                                user=args.user, layout=layout,
                                               logger_name='fasterq')

            fq_validation = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/fastqdump_validation.csv')
            fq_validation = pd.concat([fq_validation, statistics_df])
            fq_validation.to_csv('/uftp/Blood/db_calc_pipeline/tmp/fastqdump_validation.csv', index=False)

            if not 'FAILED' in statistics_df.Mark.values:
                MergeFastqFiles(df_stats=statistics_df, 
                                         sample_path=sample_path_local, 
                                         sample_id=sample_id, dataset=dataset,
                                         layout=layout, logger_name='merging')
              
                if not Path(sample_path_uftp, 'hla').exists():
                    output_optitype_path = Path('/uftp/Blood/db_calc_pipeline/tmp/optitype_validation.csv')
                    RunOptitype(sample_path=sample_path_local,
                                output_path=output_optitype_path,
                                sample_id=sample_id,
                                dataset=dataset,
                                layout=layout,
                                logger_name='optitype')

                if not Path(sample_path_uftp, 'results').exists():
                    output_mixcr_path = Path('/uftp/Blood/db_calc_pipeline/tmp/mixcr_validation.csv')
                    RunMixcr(n_start=args.n_start,
                             n_end=args.n_end,
                             sample_path=sample_path_local,
                             output_path=output_mixcr_path,
                             sample_id=sample_id,
                             dataset=dataset,
                             layout=layout,
                             logger_name='mixcr')

                if layout == 'PAIRED':
                    subprocess.run('rm {}'.format(str(Path(sample_path_local, 'sample_1.fastq'))), shell=True)
                    subprocess.run('rm {}'.format(str(Path(sample_path_local, 'sample_2.fastq'))), shell=True)
                else:
                    subprocess.run('rm {}'.format(str(Path(sample_path_local, 'sample.fastq'))), shell=True)

                subprocess.run('cp -r {} {}'.format(str(sample_path_local), str(dataset_path)),
                    shell=True, check=True)
                
                logger.info('results copied to uftp')
                subprocess.run('rm -r {}'.format(str(sample_path_local)),
                    shell=True)
            break
        if sample_id == last_SRX:
            break
    logger.info('End of calculation list')