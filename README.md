RE-Seq work flow
=====
RE-Seq is an automated workflow for analyzing RNA editing events, encompassing data quality control, alignment, acquisition of RNA editing information, multi-level statistical analysis of editing information, multi-level analysis of RNA editing events, relevant visualization analysis, and automatic report generation.

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
