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
* TEST
  * python all.py --work_dir ./01_work --out_dir ./02_Result -ref rn6 -i sample.xls -cn test -pn "Rattus norvegicus 9 RE-seq"

