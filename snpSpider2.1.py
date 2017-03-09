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

geneResult_file = open(geneResult_filename, mode='w+')
gwasResult_file = open(gwasResult_filename, mode='w+')
combinedGene_file = open(combinedGene_filename, mode='w+')
combinedGWAS_file = open(combinedGWAS_filename, mode='w+')
eqtlResult_file = open(eqtlResult_filename, mode='w+')


snpInfo_list = []
geneInfo_list = []
rsSNP_list = []
gwasTmp = []
for line in rsSNP_file.readlines():
    rsSNP_list.append(line.split()[0])
rsSNP_file.close()

print 'The rsSNPfile to be investigated: '+rsSNP_filename
print 'The genome assembly used: hg'+str(genomeAssembly)
print 'The SNPs to be investigated: '+', '.join(rsSNP_list)
print '--------------------------------------'

for rsSNP in rsSNP_list:
    print 'processing '+rsSNP+'...'
    print '--  --  --  --  --  --  --  --  --  --'
    snpTmp = snpClass(rsSNP)
    snpTmp.rseqtl = eqtlGene(rsSNP)
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
    lowerBP_snp = int(strTmp[2])
    upperBP_snp = int(strTmp[3])
    urlSNP = strTmp[0] + 'chr' + strTmp[1] + ':' + str(lowerBP_snp - 250000 + 250) + '-' + str(upperBP_snp + 250000 - 250) + strTmp[4]


    print 'opening url... '+ urlSNP
    start_ucsc = datetime.datetime.now()
    # driver = webdriver.Chrome()
    driver = webdriver.PhantomJS()
    driver.get(urlSNP)
    element = driver.find_element_by_name("hgt.hideAll")
    element.send_keys(Keys.ENTER)
    select = Select(driver.find_element_by_name('refGene'))
    select.select_by_visible_text("full")
    select = Select(driver.find_element_by_name('knownGene'))
    select.select_by_visible_text("full")
    select = Select(driver.find_element_by_name('gwasCatalog'))
    select.select_by_visible_text("full")
    select = Select(driver.find_element_by_name('snp147'))
    select.select_by_visible_text("full")
    element = driver.find_element_by_name("hgt.refresh")
    element.send_keys(Keys.ENTER)

    try:
        wait = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, '//*[@id="td_data_snp147"]/div[2]/map'))
        )
    except:
        driver.quit()
        print 'Cannot load the ucsc webpage of '+rsSNP
        sys.exit(1)

    url_Pass = driver.current_url
    rTmp = requests.get(url_Pass)
    soup = BeautifulSoup(rTmp.content, "lxml")

    # Get the band of current investigated snp
    rssnpURL = soup.find(attrs={"name":"map_data_snp147"}).find(attrs={"title":rsSNP})["href"]
    if rssnpURL[0:2] == '..':
        urlRS = 'https://genome.ucsc.edu/' + rssnpURL[3:]
    else:
        urlRS = 'https://' + rssnpURL
    requestRS = requests.get(urlRS)
    bandTmp = ' | '.join(trimStr(re.findall("<B>Bands?:</B>(.+)<BR>", requestRS.content)))
    snpTmp.band = bandTmp

    end_ucsc = datetime.datetime.now()
    print 'Time for initializing the UCSC webpage: '+str((end_ucsc-start_ucsc).seconds)+'s'

    genes = []
    recordGene_url = []
    print '--  --  --  --  --  --  --  --  --  --'
    print 'Crawling the gene data...'
    start_gene = datetime.datetime.now()
    eleGeneMap = soup.find(attrs={"name": "map_data_refGene"})
    if str(eleGeneMap) != 'None':
        for eleGene in eleGeneMap:
            if not unicode(eleGene) == '\n':
                titleGene = eleGene.attrs.get('title')
                urlGeneTmp = eleGene.attrs.get('href')
                if urlGeneTmp not in recordGene_url:
                    recordGene_url.append(urlGeneTmp)
                    if urlGeneTmp[0:2] == '..':
                        urlGene = 'https://genome.ucsc.edu/' + urlGeneTmp[3:]
                    else:
                        urlGene = 'https://' + urlGeneTmp
                    genes.append([titleGene, urlGene])

    geneInfo = []
    if len(genes) == 0:
        print 'No genes found within 500kbp'
    else:
        for ele in genes: # ele[0]: gene Title; ele[1]: gene url
            if len(re.findall('Next Exon',ele[0]))==0 and len(re.findall('Prev Exon',ele[0]))==0:
                requestGene = requests.get(ele[1])

                name = ele[0]
                description = ' | '.join(trimStr(re.findall("<B>Description:</B>(.+)<BR>", requestGene.content)))
                status = ' | '.join(trimStr(re.findall("<B>Status: </B>(.+)<BR>", requestGene.content)))
                refSeq = ' | '.join(trimStr(re.findall("<B>Status: </B>(.+)<BR>", requestGene.content)))
                band = ' | '.join(trimStr(re.findall("<B>Bands?:</B>(.+)<BR>", requestGene.content)))
                position = ' | '.join(trimStr(re.findall("<B>Position:</B> <A HREF=.+>(.+)</A><BR>", requestGene.content)))
                geneInfo.append({'name':name, 'description':description, 'status':status, 'refSeq':refSeq,
                                 'band':band, 'position':position})
                if (isInGene(position, lowerBP_snp, upperBP_snp) and (ele[0] not in rsSNP_Genes)):
                    rsSNP_Genes.append(name)

    snpTmp.geneInfo = geneInfo
    snpTmp.rsgene = rsSNP_Genes

    end_gene = datetime.datetime.now()
    print 'Time used: '+str((end_gene-start_gene).seconds)+'s'
    print '--  --  --  --  --  --  --  --  --  --'
    print 'Crawling the GWAS data...'
    start_gwas = datetime.datetime.now()
    gwasSNP_local = []
    recordGWAS_url_local = []
    eleGWASMap = soup.find(attrs={"name": "map_data_gwasCatalog"})
    if str(eleGWASMap) != 'None':
        for eleGWAS in eleGWASMap:
            if not unicode(eleGWAS) == '\n':
                titleGWAS = eleGWAS.attrs.get('title')
                urlGWASTmp = eleGWAS.attrs.get('href')
                if urlGWASTmp not in recordGWAS_url_local:
                    recordGWAS_url_local.append(urlGWASTmp)
                    if urlGWASTmp[0:2] == '..':
                        urlGWAS = 'https://genome.ucsc.edu/' + urlGWASTmp[3:]
                    else:
                        urlGWAS = 'https://' + urlGWASTmp

                    gwasSNP_local.append([titleGWAS, urlGWAS])

    end_gwas = datetime.datetime.now()
    print 'Time used: '+str((end_gwas-start_gwas).seconds)+'s'
    print '--  --  --  --  --  --  --  --  --  --'

    gwasInfo = []
    if len(gwasSNP_local)==0:
        print 'No GWAS hits found within 500kbp'
    else:
        for ele in gwasSNP_local:
            if ele[1] not in recordGWAS_url_local:
                recordGWAS_url_local.append(ele[1])
                print 'Find GWAS hits within 500kbp: '+ele[0]
                if ele[0] not in [x['name'] for x in gwasTmp]:
                    requestGWAS = requests.get(ele[1])
                    name = ele[0]
                    band = ' | '.join(trimStr(re.findall("<B>Bands?:</B>(.+)<BR>", requestGWAS.content)))
                    position = ' | '.join(trimStr(re.findall("<B>Position:</B> <A HREF=.+>(.+)</A><BR>", requestGWAS.content)))
                    disease = ' | '.join(trimStr(re.findall("<B>Disease or trait:</B>(.+)<BR>", requestGWAS.content)))
                    sample = ' | '.join(trimStr(re.findall("<B>Initial sample size:</B>(.+)<BR>", requestGWAS.content)))
                    replication = ' | '.join(trimStr(re.findall("<B>Replication sample size:</B>(.+)<BR>", requestGWAS.content)))
                    gene = ' | '.join(trimStr(re.findall("<B>Reported gene.+:</B>(.+)<BR>", requestGWAS.content)))
                    riskAllele = ' | '.join(trimStr(re.findall("<B>Reported gene.+:</B>(.+)<BR>", requestGWAS.content)))
                    pValue = ' | '.join(trimStr(re.findall("<B>p-Value:</B>(.+)<BR>", requestGWAS.content)))
                    oddsRatio =  ' | '.join(trimStr(re.findall("<B>Odds Ratio or beta:</B>(.+)<BR>", requestGWAS.content)))
                    CI = ' | '.join(trimStr(re.findall("<B>95% confidence interval:</B>(.+)<BR>", requestGWAS.content)))
                    publication = ' | '.join(trimStr(re.findall("<B>Publication:</B>(.+)<BR>", requestGWAS.content)))
                    eqtlTmp = eqtlGene(name)
                    eleGWASInfo = {'name':name, 'band':band, 'position':position, 'disease':disease, 'sample':sample,
                                     'replication':replication, 'gene':gene, 'riskAllele':riskAllele, 'pValue':pValue,
                                     'oddsRatio':oddsRatio, 'CI':CI, 'publication':publication, 'eqtl': eqtlTmp}
                    gwasInfo.append(eleGWASInfo)
                    gwasTmp.append({'name':name, 'url':ele[1], 'gwasInfo':eleGWASInfo})
                else:
                    gwasInfo.append(filter(lambda etmp: etmp['name']==ele[0], gwasTmp)[0]['gwasInfo'])
                    print 'Information added'
                    print '--  --  --  --  --  --  --  --  --  --'
    snpTmp.gwasInfo = gwasInfo
    snpInfo_list.append(snpTmp)
    driver.close()
    end = datetime.datetime.now()
    print 'Time for '+rsSNP+': '+str((end-start).seconds)+'s'
    print '--------------------------------------'
    allInfo = combineInfo(snpInfo_list)

# # Write output results
# Write headers
geneResult_file.write('*******************************************************\n'+
                        'The rsSNPfile investigated: '+rsSNP_filename+'\n'
                        +'The genome assembly used: hg'+str(genomeAssembly)+'\n'
                        +'The SNPs investigated: '+', '.join(rsSNP_list)+'\n\n')
gwasResult_file.write('*******************************************************\n'+
                        'The rsSNPfile investigated: '+rsSNP_filename+'\n'
                        +'The genome assembly used: hg'+str(genomeAssembly)+'\n'
                        +'The SNPs investigated: '+', '.join(rsSNP_list)+'\n\n')
combinedGWAS_file.write('*******************************************************\n'+
                        'The rsSNPfile investigated: '+rsSNP_filename+'\n'
                        +'The genome assembly used: hg'+str(genomeAssembly)+'\n'
                        +'The SNPs investigated: '+', '.join(rsSNP_list)+'\n\n')
combinedGene_file.write('*******************************************************\n'+
                        'The rsSNPfile investigated: '+rsSNP_filename+'\n'
                        +'The genome assembly used: hg'+str(genomeAssembly)+'\n'
                        +'The SNPs investigated: '+', '.join(rsSNP_list)+'\n\n')
eqtlResult_file.write('*******************************************************\n'+
                        'The rsSNPfile investigated: '+rsSNP_filename+'\n'
                        +'The genome assembly used: hg'+str(genomeAssembly)+'\n'
                        +'The SNPs investigated: '+', '.join(rsSNP_list)+'\n\n')
combinedGWAS_file.write('*******************************************************\n'
                            +'The Bands Involved: '+', '.join([x['band'] for x in allInfo])+'\n\n')
combinedGene_file.write('*******************************************************\n'
                            +'The Bands Involved: '+', '.join([x['band'] for x in allInfo])+'\n\n')

# Write the gene and gwas results
for eleRS in snpInfo_list:
    # Write the gene results
    geneResult_file.write('*******************************************************\n')
    geneResult_file.write('The investigated SNP: '+eleRS.name+'\n')
    geneResult_file.write('Overlapped gene(s): ')
    if len(set(eleRS.rsgene))==0:
        geneResult_file.write('No reported\n')
    else:
        geneResult_file.write(', '.join(set(eleRS.rsgene))+'\n')
    geneResult_file.write('The eQTL(s): ')
    geneResult_file.write(', '.join(set([x['gene'] for x in eleRS.rseqtl]))+'\n\n')
    geneResult_file.write('Gene(s) within the 500 kbp region: \n')
    for eleGene in eleRS.geneInfo:
        geneResult_file.write('\n<' + eleGene['name'] + '>:\n')
        geneResult_file.write('Description: ' + eleGene['description'] + '\n')
        geneResult_file.write('Status: ' + eleGene['status'] + '\n')
        geneResult_file.write('RefSeq: ' + eleGene['refSeq'] + '\n')
        geneResult_file.write('Band: ' + eleGene['band'] + '\n')
        geneResult_file.write('Position: ' + eleGene['position'] + '\n')

    # Write the gwas results
    gwasResult_file.write('*******************************************************\n')
    gwasResult_file.write('The investigated SNP: '+eleRS.name+'\n')
    gwasResult_file.write('Overlapped gene(s): ')
    if len(set(eleRS.rsgene))==0:
        gwasResult_file.write('No reported\n')
    else:
        gwasResult_file.write(', '.join(set(eleRS.rsgene))+'\n')
    gwasResult_file.write('The eQTL(s): ')
    gwasResult_file.write(', '.join(set([x['gene'] for x in eleRS.rseqtl]))+'\n')
    for eleGWAS in eleRS.gwasInfo:
        gwasResult_file.write('\n<' +eleGWAS['name'] + '>:\n')
        gwasResult_file.write('Band: ' + eleGWAS['band'] + '\n')
        gwasResult_file.write('Position: ' + eleGWAS['position'] + '\n')
        gwasResult_file.write('Disease or Trait: ' + eleGWAS['disease'] + '\n')
        gwasResult_file.write('Initial Sample Size: ' + eleGWAS['sample'] + '\n')
        gwasResult_file.write('Replication Sample Size: ' + eleGWAS['replication'] + '\n')
        gwasResult_file.write('Reported Gene(s): ' + eleGWAS['gene'] + '\n')
        gwasResult_file.write('Risk Allele Frequency: ' + eleGWAS['riskAllele'] + '\n')
        gwasResult_file.write('p-Value: ' + eleGWAS['pValue'] + '\n')
        gwasResult_file.write('Odds Ratio or Beta: ' + eleGWAS['oddsRatio'] + '\n')
        gwasResult_file.write('95% confidence interval: ' + eleGWAS['CI'] + '\n')
        gwasResult_file.write('Publication: ' + eleGWAS['publication'] + '\n')

# Write the combined gene and gwas results based on bands
for eleBand in allInfo:
    combinedGWAS_file.write('*******************************************************\n')
    combinedGWAS_file.write('The investigated Band: '+eleBand['band']+'\n')
    combinedGWAS_file.write('\n'+'The rsSNPs included in this band: \n')
    i = 0
    for eleSNP in eleBand['snp']:
        i = i+1
        combinedGWAS_file.write('('+str(i)+') '+eleSNP+': Gene(s):')
        geneTmp = filter(lambda etmp: etmp.name == eleSNP, snpInfo_list)[0].rsgene
        if len(geneTmp) == 0:
            combinedGWAS_file.write('No reported')
        else:
            combinedGWAS_file.write(', '.join(geneTmp))
        eqtlTmp = filter(lambda etmp: etmp.name == eleSNP, snpInfo_list)[0].rseqtl
        combinedGWAS_file.write('; eQTL(s):'+', '.join(set(x['gene'] for x in eqtlTmp))+'\n')
    combinedGWAS_file.write('\nGWAS hit(s) within the 500 kbp region:\n\n')
    for eleGWAS in eleBand['gwas']:
        combinedGWAS_file.write('<'+eleGWAS['name']+'>:\n')
        combinedGWAS_file.write('Band: '+eleGWAS['band']+'\n')
        combinedGWAS_file.write(('Disease or Trait: ' + eleGWAS['disease'] + '\n'))
        combinedGWAS_file.write(('Initial Sample Size: ' + eleGWAS['sample'] + '\n'))
        combinedGWAS_file.write(('Replication Sample Size: ' + eleGWAS['replication'] + '\n'))
        combinedGWAS_file.write(('Reported Gene(s): ' + eleGWAS['gene'] + '\n'))
        combinedGWAS_file.write(('Risk Allele Frequency: ' + eleGWAS['riskAllele'] + '\n'))
        combinedGWAS_file.write(('p-Value: ' + eleGWAS['pValue'] + '\n'))
        combinedGWAS_file.write(('Odds Ratio or Beta: ' + eleGWAS['oddsRatio'] + '\n'))
        combinedGWAS_file.write(('95% confidence interval: ' + eleGWAS['CI'] + '\n'))
        combinedGWAS_file.write(('Publication: ' + eleGWAS['publication'] + '\n'))
        combinedGWAS_file.write(('eQTL(s): ' + ', '.join(set([x['gene'] for x in eleGWAS['eqtl']]))+ '\n\n'))

    # Write the combined GWAS results
    combinedGene_file.write('*******************************************************\n')
    combinedGene_file.write('The investigated Band: '+eleBand['band']+'\n')
    combinedGene_file.write('\n'+'The rsSNPs included in this band: \n')
    i = 0
    for eleSNP in eleBand['snp']:
        i = i+1
        combinedGene_file.write('('+str(i)+')' +eleSNP+': Gene(s):')
        geneTmp = filter(lambda etmp: etmp.name == eleSNP, snpInfo_list)[0].rsgene
        if len(geneTmp) == 0:
            combinedGene_file.write('No reported')
        else:
            combinedGene_file.write(', '.join(geneTmp))
        eqtlTmp = filter(lambda etmp: etmp.name == eleSNP, snpInfo_list)[0].rseqtl
        combinedGene_file.write('; eQTL(s):'+', '.join(set(x['gene'] for x in eqtlTmp))+'\n')

    combinedGene_file.write('\nGene(s) within the 500 kbp region: '+', '.join(set([x['name'] for x in eleBand['gene']]))
                            +'\n\nDetailed Gene Information:\n\n')
    for eleGene in eleBand['gene']:
        combinedGene_file.write('<'+eleGene['name']+'>:\n')
        combinedGene_file.write('Description: '+eleGene['description']+'\n')
        combinedGene_file.write(('Status: ' + eleGene['status'] + '\n'))
        combinedGene_file.write(('RefSeq: ' + eleGene['refSeq'] + '\n'))
        combinedGene_file.write(('Band: ' + eleGene['band'] + '\n'))
        combinedGene_file.write(('Position: ' + eleGene['position'] + '\n\n'))

# Write the information of eQTLs
eqtlResult_file.write('*******************************************************\n')
eqtlResult_file.write('The eQTL(s) of the investigated SNPs')
for ele in snpInfo_list:
    eqtlResult_file.write('\n\n<'+ele.name+'>')
    if ele.rseqtl[0]['exist'] == False:
        eqtlResult_file.write(': '+ele.rseqtl[0]['gene']+'\n')
    else:
        i = 0
        for eleEQTL in ele.rseqtl:
            i = i+1
            eqtlResult_file.write('\n('+str(i)+') '+eleEQTL['gene']+', p Value: '+eleEQTL['pValue']
                                  +', effect size: '+eleEQTL['effectSize']+', tissue: '+eleEQTL['tissue'])

eqtlResult_file.write('\n\n*******************************************************\n')
eqtlResult_file.write('The eQTL(s) of the GWAS hits')
for ele in gwasTmp:
    eqtlResult_file.write('\n\n<'+ele['name']+'>')
    if ele['gwasInfo']['eqtl'][0]['exist'] == False:
        eqtlResult_file.write(': '+ele['gwasInfo']['eqtl'][0]['gene']+'\n')
    else:
        i = 0
        for eleEQTL in ele['gwasInfo']['eqtl']:
            i = i+1
            eqtlResult_file.write('\n('+str(i)+') '+eleEQTL['gene']+', p Value: '+eleEQTL['pValue']
                                  +', effect size: '+eleEQTL['effectSize']+', tissue: '+eleEQTL['tissue'])

geneResult_file.close()
gwasResult_file.close()
combinedGWAS_file.close()
combinedGene_file.close()
eqtlResult_file.close()

print 'The result file '+geneResult_filename+' has been output'
print 'The result file '+gwasResult_filename+' has been output'
print 'The result file '+combinedGene_filename+' has been output'
print 'The result file '+combinedGWAS_filename+' has been output'
print 'The result file '+eqtlResult_filename+' has been output'