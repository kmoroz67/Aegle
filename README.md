

## Введение
Данный пайплайн представляет собой инструмент для определения клонотипов T и B клеточных рецепторов, 
а также основных типов аллелей HLA в образце по данным РНК-секвенирования.

Для этого в пайплайн внедрены такие тулы как:

- [fasterq-dump](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump) - инструмент для конвертации  SRA данных в .fastq формат; 
является более быстрой версией инструмента fastq-dump

- [MiXCR](https://github.com/milaboratory/mixcr) - универсальный инструмент для быстрого и точного анализа T- и B клеточного репертуара на основе данных NGS.

- [Optitype](https://github.com/FRED-2/OptiType) - алгоритм для HLA генотипирования на основе данных NGS. Позволяет определять аллели HLA 1-го класса.

Подробное описание работы пайплайна доступно по [ссылке](https://bostongene.atlassian.net/wiki/spaces/BLOOD/pages/2569896111/Aesculapius+pipeline)

## Структура пайплайна:  
Пайплайн разбит на несколько файлов:

- functions.py, в котором собраны все используемые в пайплайне функции
- geo_downloader.py, в котором находится тело пайплайна по скачиванию и процессингу данных из базы данных GEO (Gene Expression Omnibus)
- tgca_downloader.py, в котором находится тело пайплайна по скачиванию и процессингу данных из базы данных GDC (Genomic Data Commons)
- main.py, из которого в зависимости от передаваемого параметра вызывается либо geo_downloader.py, либо tcga_downloader.py

## Предварительные настройки

**2. Настроить работу fasterq-dump**  

Скачать [архив](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-centos_linux64.tar.gz) и распаковать на локальной машине в папку ncbi.  
	
	$ tar -xvzf sratoolkit.2.10.9-centos_linux64-cloud.tar.gz -C ./ncbi
	
    $ docker pull fred2/optitype
	
**4. Установить инструмент samtools**

	$ apt update
	$ apt install samtools

## Запуск вычислений
#### 1. Склонировать репозиторий на локальную машину

    $ git clone https://github.com/kmoroz67/Aegle.git

#### 2. Зайти в папку проекта

    $ cd aesculapius/

#### 3. Запустить пайплайн  

Пример запуска скрипта из консоли:

    $ python3.7 main.py -u isemenov --dbase tcga --n_start 0 --n_end 7

Получить краткую справку об используемых параметрах можно c помощью следующей команды:

    $ python3.7 main.py -h

## Консольный интерфейс инструмента
Данный инструмент представляет из себя консольную утиллиту, которая поддерживает следующие параметры:

- -u / --user (опицональное) название юзера, в домашней папке которого будет создаваться папка mixcr_calc для хранения промежуточных результатов. 
	По умолчанию указывается название юзера, от имени которого запускается пайплайн.
	
- --dbase - название базы данных, из которой скачиваются и процессируются данные: geo или tcga
	
- --n_start - номер первого ядра из диапазона ядер, которые будут задействованы под расчеты MiXCR

- --n_end - номер последнего ядра из диапазона ядер, которые будут задействованы под расчеты MiXCR

	Так же доступен флаг -h / --help, отображающий краткую справку об использовании инструмента

## Хранение данных

Изначально данные считаются на локальной машине в папке calc_mixcr, затем переносятся на /uftp

Абсолютный путь: **/uftp/Blood/db_calc_pipeline**

Скаченные данные имеют следующую структуру:  dataset --> sample --> hla, results  

- hla - папка с результатами работы Optitype  

- results - папка с результатами работы MiXCR

Также имеется папка **tmp**, в которой хранятся следующие файлы:  

Для данных с GEO:

- **ann_calc_slice.csv** - срез таблицы Annotation calculated   

	Run - название рана  
	Sample  - название образца  
	Dataset - название датасета  
	Layout  - тип рана (single или paired)  
	n_repeats - количество попыток запустить fasterq-dump и скачать .fastq файлы  
	
- **download_process_samples.csv** - список образцов, которые попали в очередь на обработку пайплайном

	Sample - название образца  

- **fastqdump_validation.csv**

	Run - название рана   
	Sample - название образца  
	Dataset - название датасета  
	Mark - оценка того, удалось ли скачать .fastq файлы или нет  
	N_repeats - количество попыток запустить fasterq-dump и скачать .fastq файлы  
	
- **optitype_validation.csv**  

 	Sample - название образца  
	Dataset - название датасета  
	Label - тип рана (single или paired)  
	Mark - оценка того, удачно ди был запущен Optitype

- **mixcr_validation.csv**

	Sample - название образца  
	Dataset - название датасета  
	Layout - тип рана (single или paired)  
	Mark  - оценка того, удачно ли был запущен MiXCR
	
- **pipeline.log**
	
	Содержит в себе логи по работе файла geo_downloader.py
	
Для данных с TCGA аналогичные файлы, только с приставкой tcga_...
