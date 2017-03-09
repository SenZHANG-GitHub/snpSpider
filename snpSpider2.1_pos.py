from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.common.by import By
from selenium.webdriver.support import expected_conditions as EC
from bs4 import BeautifulSoup
from collections import Counter
import datetime, requests, re, sys, os, errno

class snpClass:
    def __init__(self, name):
        self.name = name
        self.bp = ''
        self.rsgene = []
        self.rseqtl = []
        self.geneInfo = []
        self.gwasInfo = []
        self.band = ''

class logger(object):
  def __init__(self, file):
    self.terminal = sys.stdout
    dir = os.path.dirname(file)
    if dir=='': dir='.'
    if not os.path.exists(dir):
      try:
        os.makedirs(os.path.dirname(file))
      except OSError as exc: # Guard against race condition
        if exc.errno != errno.EEXIST: raise
    self.log = open(file, "w")

  def write(self, message):
    self.terminal.write(message)
    self.log.write(message)

  def flush(self):
    #this flush method is needed for python 3 compatibility.
    #this handles the flush command by doing nothing.
    #you might want to specify some extra behavior here.
    pass

  def __del__(self):
    self.log.close()

def trimStr(t):
    tNew = []
    for ele in t:
        ele = re.sub("<em>", '', ele)
        ele = re.sub("</em>", '', ele)
        ele = re.sub("</A>", '', ele)
        ele = re.sub("<A HREF.+>", '', ele)
        tNew.append(ele)
    return tNew

def isInGene(t,lowerBP,upperBP):
    if len(t) == 0:
        return False
    else:
        tN = re.findall("chr.+:(\d+)-(\d+)",''.join(t))
        lowerGene = int(tN[0][0])
        upperGene = int(tN[0][1])
        if not (upperBP < lowerGene or upperGene < lowerBP):
            return True
        else:
            return False

def combineInfo(snpInfo):
    bandList = list(set([x.band for x in snpInfo]))
    gwasList = [[] for i in bandList]
    snpOri = [[] for i in bandList]
    geneList = [[] for i in bandList]
    for ele in snpInfo:
        snpOri[bandList.index(ele.band)].append(ele.name)
        for gsnp in ele.gwasInfo:
            if gsnp['name'] not in [x['name'] for x in gwasList[bandList.index(ele.band)]]:
                gwasList[bandList.index(ele.band)].append(gsnp)
        for ggene in ele.geneInfo:
            if ggene not in geneList[bandList.index(ele.band)]:
                geneList[bandList.index(ele.band)].append(ggene)
    allInfo = []
    for ele in bandList:
        allInfo.append({'band':ele, 'snp': snpOri[bandList.index(ele)], 'gwas':gwasList[bandList.index(ele)],
                       'gene':geneList[bandList.index(ele)]})
    return allInfo

def eqtlGene(rs):
    eQTL_list = []
    start_eqtl = datetime.datetime.now()

    urlEQTL = 'http://www.gtexportal.org/home/eqtls/bySnp?snpId=' + rs+'&tissueName=All'
    print 'opening url... '+ urlEQTL
    driverEQTL = webdriver.PhantomJS()
    driverEQTL.get(urlEQTL)

    try:
        wait = WebDriverWait(driverEQTL, 10).until(
            EC.presence_of_element_located((By.XPATH, '//*[@id="gridTable"]/tbody'))
        )
    except:
        driverEQTL.quit()
        print 'Cannot load the eqtl webpage of '+rs
        sys.exit(1)

    sourceTmp = driverEQTL.page_source
    sourceTrim = re.sub("\n", " ", sourceTmp)
    reRule = "<div id=\"singleEqtl\".+<div id=\"tableDiv\">.+<table id=\"gridTable\".+<tbody>+(.+)</tbody>"
    geneUnformated = re.findall(reRule, sourceTrim.split("eGeneTitle")[0])[0].split("</tr>")
    geneUnformated.pop()
    for ele in geneUnformated:
        if len(re.findall('No significant eQTLs',geneUnformated[0])) > 0:
            exist = False
            gene = 'No significant eQTLs'
            eQTL_list.append({'exist': exist, 'gene': gene})
        else:
            exist = True
            eleRows = ele.split("</td>")
            gene = re.findall('<a href.+>(.+)</a>', eleRows[1])[0]
            pValue = re.findall('<td>(.+)', eleRows[3])[0]
            effectSize = re.findall('<td>(.+)', eleRows[4])[0]
            tissue = re.findall('<a href.+>(.+)</a>', eleRows[5])[0]
            eQTL_list.append({'exist': exist, 'gene': gene, 'pValue': pValue, 'effectSize': effectSize, 'tissue': tissue})
    end_eqtl = datetime.datetime.now()
    driverEQTL.close()
    print 'Time for investigating the eQTLs of '+rs+': '+str((end_eqtl-start_eqtl).seconds)+'s'
    print '--  --  --  --  --  --  --  --  --  --'
    return eQTL_list


print '==============snpSpider==============='

helpText = '''
  =================HELP=================
  --snplist: SNP list file (Default: 'rsSNPList.txt')
  --assembly: genome assembly (Default: 38, Corresponding to hg38) (Options: 38, 19, 18, 17, 16)
  --out: outfilenames followed by default names (Default: results_Gene.txt, results_GWAS.txt,results_combinedGWAS.txt, results_combinedGene.txt, results_EQTL.txt)
  --log: log file (Default: log.txt)
  --h: Help
  ======================================'''

arg={'--snplist':None, '--assembly':None, '--out': None, '--log': 'log.txt', '--h': None}

if len(sys.argv)%2==1:
  for i in range(len(sys.argv)/2):
    arg[sys.argv[2*i+1]]=sys.argv[2*i+2]
else:
  if sys.argv[1]!='--h': print 'The number of arguments is wrong!'
  print helpText
  sys.exit(1)

sys.stdout = logger(arg['--log'])
  
if arg['--snplist'] != None:
    rsSNP_filename = arg['--snplist']
else:
    rsSNP_filename = 'rsSNPList.txt'
	
if arg['--assembly'] != None:
    genomeAssembly = int(arg['--assembly'])
else:
    genomeAssembly = 38

if arg['--out'] != None:
    geneResult_filename = arg['--out']+'_results_Gene.txt'
    gwasResult_filename = arg['--out']+'_results_GWAS.txt'
    combinedGene_filename = arg['--out']+'_results_combinedGene.txt'
    combinedGWAS_filename = arg['--out']+'_results_combinedGWAS.txt'
    eqtlResult_filename = arg['--out']+'_results_EQTL.txt'
else:
    geneResult_filename = 'results_Gene.txt'
    gwasResult_filename = 'results_GWAS.txt'
    combinedGene_filename = 'results_combinedGene.txt'
    combinedGWAS_filename = 'results_combinedGWAS.txt'
    eqtlResult_filename = 'results_EQTL.txt'
	
if arg['--h'] != None:
    print helpText

try:
    rsSNP_file = open(rsSNP_filename,mode='r')
except IOError:
    print 'Cound not read snp list file: '+rsSNP_filename
    sys.exit(1)


snpInfo_list = []
geneInfo_list = []
rsSNP_list = []
rsSNPbp_list = []
lowbp_list = []
upbp_list = []
gwasTmp = []
for line in rsSNP_file.readlines():
    rsSNP_list.append(line.split()[0])
rsSNP_file.close()

print 'The rsSNPfile to be investigated: '+rsSNP_filename
print 'The genome assembly used: hg'+str(genomeAssembly)
print 'The SNPs to be investigated: '+', '.join(rsSNP_list)
print '--------------------------------------'

fout = open("overlap_snpName_out.pos", mode  ="w")

for i in range(len(rsSNP_list)):
    rsSNP = rsSNP_list[i]
    print 'processing '+rsSNP+'...'
    print '--  --  --  --  --  --  --  --  --  --'
    snpTmp = snpClass(rsSNP)
    rsSNP_Genes = []
    start = datetime.datetime.now()

    # # Deal with the web url, 500kbp centered on certain SNP with certain genome assembly e.g. hg38
    # # Using rs ID to find out the corresponding ucsc url
    # urlSNP = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg'+str(genomeAssembly)+'&position=chr'+str(chr)+'%3A'+str(lowerBP)+'-'+str(upperBP)
    # Using soup.find will give out the top/first result

    urlSNPTmp = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg' + str(genomeAssembly) + '&position=' + rsSNP
    rTmp = requests.get(urlSNPTmp)
    soup = BeautifulSoup(rTmp.content, 'lxml')
    eleTmp = soup.find('div', id='hgFindResults')

    urlSNPInit = 'https://genome.ucsc.edu/cgi-bin/' + eleTmp.find('a').attrs.get('href')
    strTmp = re.split(r"chr(\d+):(\d+)-(\d+)", urlSNPInit)
    chrtmp = strTmp[1]
    lowerBP_snp = int(strTmp[2])
    upperBP_snp = int(strTmp[3])
    snpPos = str((lowerBP_snp+upperBP_snp)/2)
    snpTmp.bp = snpPos
    rsSNPbp_list.append(snpPos)
    lowbp_list.append(str(lowerBP_snp))
    upbp_list.append(str(upperBP_snp))

    # need to use str.format() for alignment
    fout.write("{:15}: \tchr {:5}, \t {:10} - {:10} \t (center: {:10})\n".format(rsSNP_list[i], chrtmp, lowbp_list[i], upbp_list[i], rsSNPbp_list[i]))
    
fout.close()


