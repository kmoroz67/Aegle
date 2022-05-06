import pandas as pd
import numpy as np
import os
from pathlib import Path
import subprocess
from itertools import islice
import shutil
import logging
import time


def SetLogger(logger_name, file_name=None):
    """
    Create custom logger and set its configuration
    
    :param logger_name: name of created logger
    """
    #Create a custom logger
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)
    
    if logger_name != 'argparse':
      
        #Create handlers
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler('/uftp/Blood/db_calc_pipeline/tmp/{file_name}.log'.format(file_name=file_name))
        c_handler.setLevel(logging.INFO)
        f_handler.setLevel(logging.WARNING)

        #Create formatters and add them to handlers
        c_format = logging.Formatter('%(levelname)s: %(message)s')
        f_format = logging.Formatter('%(asctime)s::%(name)s::%(levelname)s::%(message)s',  datefmt='%Y-%m-%d %H:%M:%S')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)

        #Add handlers to the logger
        logger.addHandler(c_handler)
        logger.addHandler(f_handler)
    
    else:
        c_handler = logging.StreamHandler()  
        c_handler.setLevel(logging.INFO)
        c_format = logging.Formatter('%(message)s')
        c_handler.setFormatter(c_format)
        logger.addHandler(c_handler)
        
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

def DownloadTCGA(run, uuid, sample_id, dataset, run_path, df, df_stats, user, logger_name='tcga'):
    """
    Download .bam files from gdc and convert them to .fastq format

    :param uuid: universally unique identifier for downloading .bam (alternative entry for run)
    :param run: name of run (TCGA-... from GDC) in calculated dataset
    :param sample_id: name of sample (TCGA-... from GDC) in calculated dataset
    :param dataset: name of dataset 
    :param run_path: path to TCGA folder (where the data will be downloaded)
    :param df: dataframe of calculation list (list of datasets to download)
    :param df_stats: cumulative dataframe with statistics for one sample
    :param user: user name to get access to fasterq-dump tool
    :param logger_name: name of created logger (by default logger_name = "TCGA")
    :return: df_stats - cumulative dataframe with statistics for one sample
    """
     
    logger = logging.getLogger(logger_name)
    
    logger.info('start download .bam file for uuid {uuid}'.format(uuid=uuid))
    subprocess.run('/uftp/shared/Software/TCGA/gdc-client download -t \
                               /uftp/shared/Datasets/TCGA/token.txt {uuid}'.format(uuid=uuid), shell = True)

    subprocess.run('cp {run_path_old}/*.bam {run_path_new}'.format(run_path_old=uuid, run_path_new=run_path), 
                   shell = True)
    subprocess.run('rm -r {run_path_old}'.format(run_path_old=uuid), shell = True)
    
    for item in run_path.iterdir():
        if item.is_file() and str(item).endswith('.bam'):
            bam_sorted_path = Path(run_path, item.stem + '_sorted.bam')
            logger.info('start sorting for {item}'.format(item=str(item.name)))
            subprocess.run('samtools sort -n {item} -o {bam_sorted_path}'.format(item=str(item), 
                                                bam_sorted_path=str(bam_sorted_path)), shell = True)
            logger.info('done sorting for {item}'.format(item=str(item.name)))
            logger.info('file {bam_sorted_path} is successfully created'.format(bam_sorted_path=str(bam_sorted_path.name)))
            
            logger.info('start checking layout...')
            is_paired = int(subprocess.check_output(['samtools', 'view', '-c', '-f', '1', str(item)]).split()[0])
            
            if is_paired:
                logger.info('PAIRED layout')
                layout = 'PAIRED'
                logger.info('start bam to fastq convertation for paired-end')
                subprocess.run('samtools fastq -1 {paired_1} -2 {paired_2} -0 /dev/null -s/dev/null {bam_sorted_path}'.format(
                            paired_1=str(Path(run_path, 'paired_1.fastq')), 
                            paired_2=str(Path(run_path, 'paired_2.fastq')),
                            bam_sorted_path=str(bam_sorted_path)), shell = True)
                logger.info('done bam to fastq convertation for paired-end')
            else:
                logger.info('SINGLE layout')
                layout = 'SINGLE' 
                logger.info('start bam to fastq convertation for single-end')
                subprocess.run('samtools fastq {bam_sorted_path} > {single}'.format(bam_sorted_path=str(bam_sorted_path),
                               single=str(Path(run_path, 'single.fastq'))), shell = True)
                logger.info('done bam to fastq convertation for single-end')
   
                   
    df_stats = pd.concat([df_stats,
                          pd.DataFrame(np.array([run, sample_id, dataset]).reshape(1, -1),
                                       columns=['Run', 'Sample', 'Dataset'])])
    logger.info('done downloading and convertation for {run}'.format(run=run))
    return (layout, df_stats)


def MergeFastqFiles(df_stats, sample_path, sample_id, dataset, layout, logger_name='merging'):
    """
    Merge all .fastq files (all runs) located inside one sample and make its validation. 
    Paired runs are merged into two separate files named sample_1.fastq and sample_2.fastq respectively.
    Single runs are merged into file named sample.fastq.

    :param df_stats: cumulative dataframe with statistics for one sample
    :param sample_path: path to sample folder (SRX for GEO, TCGA for GDC) with downloaded data
    :param sample_id: name of sample (SRX from GEO, TCGA from GDC) in calculated dataset
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
