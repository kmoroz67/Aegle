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

def CalculateTCGA(user, n_start, n_end):
    
    #create loggers
    loggers_used = ['main', 'tcga', 'merging', 'optitype', 'mixcr']
    for logger in loggers_used:
        functions.SetLogger(logger_name=logger, file_name='tcga_pipeline')
    
    #start the main cycle of list check
    logger = logging.getLogger('main')
    parent_path = Path('/uftp/Blood/db_calc_pipeline/')
    local_path = Path('/home/{user}/mixcr_calc/'.format(user=user))
    while True:
        df = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/tcga_annotation_slice.tsv', sep='\t')
        last_TCGA = df['sample.submitter_id'].tail(1).item()
        datasets = df.project.unique()
        np.random.shuffle(datasets)
        for dataset in datasets:
            dataset_path = parent_path / dataset
            dataset_path.mkdir(parents=True, exist_ok=True)
            for sample_id in df['sample.submitter_id'][df.project == dataset]:

                sample_path_uftp = parent_path / dataset / sample_id
                sample_path_local = local_path / dataset / sample_id

                if Path(sample_path_uftp, 'hla').exists() and Path(sample_path_uftp, 'results').exists():
                    need_to_calculate = False
                    continue
                else:  
                    df_process = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/tcga_download_process_samples.csv')
                    if sample_id in df_process['Sample'].values:
                        need_to_calculate = False
                        continue      
                    else:      
                        with open('/uftp/Blood/db_calc_pipeline/tmp/tcga_download_process_samples.csv', 'a') as f:
                            f.write('{}\n'.format(sample_id))
                        need_to_calculate = True
                        break  #####################################
                        
            if not need_to_calculate:
                continue
                                
            sample_path_local.mkdir(parents=True, exist_ok=True)
            subprocess.run('chmod 777 {}'.format(str(sample_path_local)), shell=True)

            statistics_df = pd.DataFrame(columns=['Run', 'Sample', 'Dataset'])
            for uuid in df.file_id[df['sample.submitter_id'] == sample_id]:
                run = df.cases[df.file_id == uuid].values[0]
                run_path = sample_path_local / run
                run_path.mkdir(parents=True, exist_ok=True)
                layout, statistics_df = functions.DownloadTCGA(run=run, uuid=uuid,
                                                          sample_id=sample_id,
                                                          dataset=dataset,
                                                          run_path=run_path, df=df, 
                                                          df_stats=statistics_df, 
                                                          user=user,
                                                          logger_name='tcga')

            fq_validation = pd.read_csv('/uftp/Blood/db_calc_pipeline/tmp/tcga_downloaded.csv')
            fq_validation = pd.concat([fq_validation, statistics_df])
            fq_validation.to_csv('/uftp/Blood/db_calc_pipeline/tmp/tcga_downloaded.csv', index=False)

            functions.MergeFastqFiles(df_stats=statistics_df, 
                                     sample_path=sample_path_local, 
                                     sample_id=sample_id, dataset=dataset,
                                     layout=layout, logger_name='merging')

            if not Path(sample_path_uftp, 'hla').exists():
                output_optitype_path = Path('/uftp/Blood/db_calc_pipeline/tmp/tcga_optitype_validation.csv')
                functions.RunOptitype(sample_path=sample_path_local,
                            output_path=output_optitype_path,
                            sample_id=sample_id,
                            dataset=dataset,
                            layout=layout,
                            logger_name='optitype')

            if not Path(sample_path_uftp, 'results').exists():
                output_mixcr_path = Path('/uftp/Blood/db_calc_pipeline/tmp/tcga_mixcr_validation.csv')
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
        if sample_id == last_TCGA:
            break
    logger.info('End of calculation list')