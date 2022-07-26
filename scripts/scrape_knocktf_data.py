## python script to scrape KnockTF database

from selenium import webdriver
import pandas as pd
import urllib3
import wget
import os
import requests
from bs4 import BeautifulSoup

def scrape_data():
    req = urllib3.PoolManager()

    URL = 'http://www.licpathway.net/KnockTF/download.php'
    download_URL = 'http://www.licpathway.net/KnockTF'
    res = req.request('GET',URL)

    soup = BeautifulSoup(res.data, 'html.parser')

    id_columns = ['Knock-Method',
     'TF',
     'TF Class',
     'TF Superclass',
     'Biosample Name',
     'TissueType',
     'Biosample Type',
     'Data Source',
     'Profile ID',
     'Platform',
     'Pubmed ID',
     'Download']

    fh_idx = [0,1,4,5,6,7,8]
    tr_tags = soup.findAll('tr')
    for tr_tag in tr_tags:
        if len(tr_tag) > 24:
            td_tags = tr_tag.select('td')
            if len(td_tags) > 11:
                try:
                    if td_tags[-1].contents[0] == 'txt ':
                        print(wget_handle)
                        wget_handle = os.path.join(download_URL, td_tags[-1].find(href=True).attrs['href'])
                        file_handle = '-'.join([td_tags[i].contents[0].replace(' ','_').replace('/','_') for i in fh_idx]) +'.tsv'
                        wget.download(wget_handle, out=file_handle )
                        print('   Complete!')
                except:
                    continue


    resp = requests.get(URL)
    soup = BeautifulSoup(resp.content, 'html5lib')

    fh_idx = [0,1,4,5,6,7,8]
    tr_tags = soup.findAll('tr')
    n = 0
    for tr_tag in tr_tags:
        if len(tr_tag) > 24:
            td_tags = tr_tag.select('td')
            if len(td_tags) > 11:
                n += 1
                if td_tags[-1].contents[0] == 'txt ':
                    wget_handle = os.path.join(download_URL, td_tags[-1].find(href=True).attrs['href'])
                    file_handle = '-'.join([td_tags[i].contents[0].replace(' ','_').replace('/','_').replace('-','_') for i in fh_idx]) +'.tsv'
                    print(file_handle)
                    wget.download(wget_handle, out=file_handle )
                    print('   Complete!')


def parse_and_generate_metadata():
    files = list(set([x.split('-')[-1].split('.tsv')[0] for x in os.listdir() if x.endswith('.tsv')]))
    for geo_id in files:
        try:
            gse = GEOparse.get_GEO(geo = geo_id, destdir="./test")
            meta = gse.phenotype_data
            out_file = 'metadata/{}_metadata.tsv'.format(gse.name)
            meta.to_csv(out_file,sep='\t')
        except:
            print('Unable to parse {}'.format(geo_id))


