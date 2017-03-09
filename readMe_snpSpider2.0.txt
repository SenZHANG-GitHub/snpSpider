****************************************************
<Author: Sen ZHANG>
<Version: snpSpider2.0>
<Please contact szhangat@connect.ust.hk for any bugs or functions needed for future development>
<The codes was developed under the Interpreter : Anaconda2 python 2.7.11>

****************************************************
<This tool is used for crawling snp information from ucsc genome browser and eQTL information from gtexportal>
<To speed up the crawling process, please turn to snpSpider1.0 without crawling the eQTL information>

****************************************************
<Need to install selenium, bs4, requests, datetime>
Enter the following commands in cmd 
"pip install selenium"
"pip install bs4"
"pip install requests"
"pip install datetime"

****************************************************
<Need to put the phantomjs executable file in PATH or in the current folder>
<Please turn to http://phantomjs.org/download.html for appropriate versions under Windows, Linux or Mac>

****************************************************
<The input format of the python script>
  =================HELP=================
  --snplist: SNP list file (Default: 'rsSNPList.txt')
  --assembly: genome assembly (Default: 38, Corresponding to hg38) (Options: 38, 19, 18, 17, 16)
  --out: outfilenames followed by default names (Default: geneResults.txt, gwasResults.txt, combined_gwasResults.txt, combined geneResults.txt, eqtlResults.txt)
  --log: log file (Default: log.txt)
  --h: Help
  ======================================

<Command line examples>
python snpSpider.py
python snpSpider.py --h
python snpSpider.py --snplist xx.txt --assembly 38 --out myresult —log mylog.txt

****************************************************
<The rsSNPs to be investigated should be put in the rsSNPList.txt(Default) file>
<An example: (The snps are put in different lines)>
rs1550586
rs11636035
rs4693299
rs35607436

