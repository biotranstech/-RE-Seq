RE-Seq work flow
=====
RE-Seq is an automated workflow for analyzing RNA editing events, written using the python DAGflow package, encompassing data quality control, alignment, acquisition of RNA editing information, multi-level statistical analysis of editing information, multi-level analysis of RNA editing events, relevant visualization analysis, and automatic report generation.

## Workflow file directory structure

```
├── all.py                    #Workflow Complete process
├── dagflow                   #python DAGflow package
├── report_noAS.html          #Report template
├── report_RNA_edic.py        #Automated report generation scripts
├── report_utils.py           #Report generating python modules
├── REseq                     #RE-Seq python modules
│   ├── common.py             #Public files handle python packages
│   ├── config.py             #Process software configuration file
│   ├── down_analysis.py      #RE-Seq downstream analysis module
│   ├── gene_count.py         #RNA expression analysis module
│   ├── qc_hista.py           #Quality control and comparison module
│   ├── report_utils.py       #Report generating python modules
│   └── work_bam_reseq.py     #Get the RNA editing information module
├── scripts                   #Analysis script path
├── template                  #html configuration file
│   ├── asset                 #html related js ccs file
│   └── images                #html related graphics files
```

Installation
-----

* REQUIREMENT
  * [Python](https://www.python.org/download/releases/)
  * [perl](https://www.perl.org)
  * [R](https://www.r-project.org)
  * [GATK](https://github.com/broadinstitute/gatk/releases)

Usage
-----
```python
usage: all.py [-h] -i STR [-ref FILE] [--trim INT] [--qvalue INT] [-m FILE] [--thread INT] [--concurrent INT] [--refresh INT]<br>
              [--job_type {sge,local}] [--work_dir DIR] [--out_dir DIR] [-as {yes,no}] -cn CONTRACT -pn NAME [-type TYPE]<br>
              [-type_num TYPE_NUM]<br>

version: v1.0.0

optional arguments:
  -h, --help            show this help message and exit
  -i STR, --input STR   Input the name of the sample and Second generation sequencing path.
  -ref FILE, --reference FILE
                        Input the host's reference database. Default: hg38. you can choose mm10 or rn6
  --trim INT            Set trim length, default=5
  --qvalue INT          The quality value that a base is qualified, default=20
  -m FILE, --method FILE
                        RE-seq method. Default: gatk. you can choose gatk; ; REDItools ; VarScan ; Sprint ; RED_ML. if you want choose
                        Multiple methods; you can use: gatk;REDItools;VarScan
  --thread INT          Analysis of the number of threads used, default=4
  --concurrent INT      Maximum number of jobs concurrent (default: 10)
  --refresh INT         Refresh time of log in seconds (default: 30)
  --job_type {sge,local}
                        Jobs run on [sge, local] (default: local)
  --work_dir DIR        Work directory (default: current directory)
  --out_dir DIR         Output directory (default: current directory)
  -as {yes,no}, --asan {yes,no}
                        AS analysis or not, yes or no, default no
  -cn CONTRACT, --contract CONTRACT
                        contract number
  -pn NAME, --name NAME
                        project name
  -type TYPE, --type TYPE
                        padj or pval
  -type_num TYPE_NUM, --type_num TYPE_NUM
                        padj or pval threshold value
```

## Status tracking and restart

The system displays the task progress and whether the task is running successfully. If the workflow interruption needs to be restarted, the system will automatically check the file, skip the completed task, and continue to execute from the breakpoint, improving work efficiency.
