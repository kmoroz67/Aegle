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

def CalculateGEO(user, n_start, n_end):
    
    #create loggers
    loggers_used = ['main', 'fasterq', 'merging', 'optitype', 'mixcr']
    for logger in loggers_used:
        functions.SetLogger(logger_name=logger, file_name='pipeline')
    
    #start the main cycle of list check
    logger = logging.getLogger('main')
    parent_path = Path('/uftp/Blood/db_calc_pipeline/')
    local_path = Path('/home/{user}/mixcr_calc/'.format(user=user))
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
                statistics_df = functions.DownloadFasterq(run=run, sample_id=sample_id,
                                                          dataset=dataset,
                                                          run_path=run_path, df=df, 
                                                          df_stats=statistics_df, 
                                                          user=user, layout=layout,
                                                          logger_name='fasterq')

            fq_validation = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/fastqdump_validation.csv')
            fq_validation = pd.concat([fq_validation, statistics_df])
            fq_validation.to_csv('/uftp/Blood/db_calc_pipeline/tmp/fastqdump_validation.csv', index=False)

            if not 'FAILED' in statistics_df.Mark.values:
                functions.MergeFastqFiles(df_stats=statistics_df, 
                                         sample_path=sample_path_local, 
                                         sample_id=sample_id, dataset=dataset,
                                         layout=layout, logger_name='merging')
              
                if not Path(sample_path_uftp, 'hla').exists():
                    output_optitype_path = Path('/uftp/Blood/db_calc_pipeline/tmp/optitype_validation.csv')
                    functions.RunOptitype(sample_path=sample_path_local,
                                output_path=output_optitype_path,
                                sample_id=sample_id,
                                dataset=dataset,
                                layout=layout,
                                logger_name='optitype')

                if not Path(sample_path_uftp, 'results').exists():
                    output_mixcr_path = Path('/uftp/Blood/db_calc_pipeline/tmp/mixcr_validation.csv')
                    functions.RunMixcr(n_start=n_start,
                             n_end=n_end,
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