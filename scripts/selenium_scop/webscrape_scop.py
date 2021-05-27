from selenium import webdriver 
from selenium.webdriver.common.by import By 
from selenium.webdriver.support.ui import WebDriverWait 
from selenium.webdriver.support import expected_conditions as EC 
from selenium.common.exceptions import TimeoutException
import time
from Bio import SeqIO
import re


#input_file = "/Users/rmylonas/Work/PAF/projects/pumba/data/fasta/tiny_UP000005640_9606.fasta"
input_file = "/Users/rmylonas/Work/PAF/projects/pumba/data/fasta/a.fasta"
output_path = "/Users/rmylonas/tmp/datamining_pumba/scop/raw/"

# parse fasta
fasta_sequences = SeqIO.parse(open(input_file),'fasta')


def run_scop_blast(sequence_str):
	browser.get('http://phosphatome.net/blast-scop/blast.html')
	sequence = browser.find_element_by_name("SEQUENCE")
	sequence.send_keys(sequence_str)
	sequence.submit()
	time.sleep(5)
	res = browser.find_element_by_xpath("(//body)")
	return res.text


def parse_ac(fasta_header):
	m = re.search(r'.+?\|(.+?)\|.+', fasta_header)
	return m.group(1)

def write_res(filepath, res):
	f = open(filepath, 'w')
	f.write(res)
	f.close()


# connect to browser
browser = webdriver.Firefox(executable_path='/Users/rmylonas/Applications/geckodriver')

for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    filepath = output_path + parse_ac(name) + ".txt"
    res = run_scop_blast(sequence)
    write_res(filepath, res)


# close the browser
browser.quit() 




