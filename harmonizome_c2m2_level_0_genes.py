import os,sys
import logging
from c2m2_frictionless import create_datapackage, validate_datapackage, validate_id_namespace_name_uniqueness
from c2m2_frictionless.C2M2_Level_0 import table_schema_specs_for_level_0_c2m2_encoding_of_dcc_metadata as c2m2_level_0
from dataclasses import dataclass, asdict
import hashlib
import requests
from tqdm import tqdm

def convert_harmonizome_genes_to_c2m2():
  ''' Construct a set of c2m2 objects corresponding to datasets available from Harmonizome.
  '''
  ns = c2m2_level_0.id_namespace(
    id='https://maayanlab.cloud/Harmonizome/',
    name='Harmonizome',
    description='Harmonizome is a collection of information about genes and proteins from datasets provided by various online resources distilled into attribute tables that define significant associations between genes and attributes, where attributes could be genes, proteins, cell lines, tissues, experimental perturbations, diseases, phenotypes, or drugs depending on the dataset.',
  )
  yield ns

  # Catalog single gene files from Harmonizome
  genes = open('Harmonizome_genes.txt').read().split('\n')
  gene_file_url = 'https://maayanlab.cloud/Harmonizome/api/1.0/download/associations?gene=%s' 
  for gene in tqdm(genes):
    try:
        with requests.get(gene_file_url % gene, stream=True) as response:
            md5 = hashlib.md5()
            sha256 = hashlib.sha256()
            filesize = 0
            for chunk in response.iter_content(4096):
                sha256.update(chunk)
                md5.update(chunk)
                filesize += len(chunk)
        sha256_hash = sha256.hexdigest()
        md5_hash = md5.hexdigest()
          
        f = c2m2_level_0.file(
        id_namespace=ns.id,
        id=gene,
        filename='%s_associations.tsv' % gene,
        size_in_bytes=filesize,
        persistent_id='https://maayanlab.cloud/Harmonizome/gene/%s' % gene,
        sha256=sha256_hash,
        md5=md5_hash,
        )
        yield f
    except Exception as e:
      logging.error('{}'.format(e))


if __name__ == '__main__':
  if len(sys.argv) >= 2:
    outdir = sys.argv[1]
  else:
    logging.error("OUTDIR must be specified.")
    sys.exit()
  pkg = create_datapackage('C2M2_Level_0', convert_harmonizome_genes_to_c2m2(), outdir)
  validate_datapackage(pkg)
  validate_id_namespace_name_uniqueness(pkg)