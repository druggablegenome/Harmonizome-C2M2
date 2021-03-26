import os,sys
import logging
from c2m2_frictionless import create_datapackage, validate_datapackage, validate_id_namespace_name_uniqueness
from c2m2_frictionless.C2M2_Level_1 import table_schema_specs_for_level_1_c2m2_encoding_of_dcc_metadata as c2m2_level_1
from dataclasses import dataclass, asdict
import hashlib
import requests
from tqdm import tqdm

def convert_harmonizome_datasets_to_c2m2():
  ''' Construct a set of c2m2 objects corresponding to datasets available from Harmonizome.
  '''
  ns = c2m2_level_1.id_namespace(
    id='http://www.druggablegenome.net/',
    abbreviation='IDG',
    name='IDG: Illuminating the Druggable Genome',
    description='The goal of IDG is to improve our understanding of understudied proteins from the three most commonly drug-targeted protein families: G-protein coupled receptors, ion channels, and protein kinases.',
  )
  yield ns

  project = c2m2_level_1.project(
    id_namespace=ns.id,
    local_id='Harmonizome',
    persistent_id='https://maayanlab.cloud/Harmonizome/',
    description='A collection of information about genes and proteins from 114 datasets provided by 66 online resources.'

  )
  yield project

  selected_datasets = [
    ('Achilles Cell Line Gene Essentiality Profiles', 'achilles'),
    ('Allen Brain Atlas Adult Human Brain Tissue Gene Expression Profiles', 'brainatlasadulthuman'),
    ('Allen Brain Atlas Adult Mouse Brain Tissue Gene Expression Profiles', 'brainatlasadultmouse'),
    ('Allen Brain Atlas Developing Human Brain Tissue Gene Expression Profiles by Microarray', 'brainatlasdevelopmentalhumanmicroarray'),
    ('Allen Brain Atlas Developing Human Brain Tissue Gene Expression Profiles by RNA-seq', 'brainatlasdevelopmentalhumanrnaseq'),
    ('Allen Brain Atlas Prenatal Human Brain Tissue Gene Expression Profiles', 'brainatlasprenatalhuman'),
    ('BIND Biomolecular Interactions', 'bind'),
    ('BioGPS Cell Line Gene Expression Profiles', 'biogpsnci60'),
    ('BioGPS Human Cell Type and Tissue Gene Expression Profiles', 'biogpshuman'),
    ('BioGPS Mouse Cell Type and Tissue Gene Expression Profiles', 'biogpsmouse'),
    ('BioGRID Protein-Protein Interactions', 'biogrid'),
    ('Biocarta Pathways', 'biocarta'),
    ('CCLE Cell Line Gene CNV Profiles', 'cclecnv'),
    ('CCLE Cell Line Gene Expression Profiles', 'cclemrna'),
    ('CCLE Cell Line Gene Mutation Profiles', 'cclemut'),
    ('CHEA Transcription Factor Binding Site Profiles', 'chea'),
    ('CHEA Transcription Factor Targets', 'cheappi'),
    ('CMAP Signatures of Differentially Expressed Genes for Small Molecules', 'cmap'),
    ('COMPARTMENTS Curated Protein Localization Evidence Scores', 'jensencompartmentcurated'),
    ('COMPARTMENTS Experimental Protein Localization Evidence Scores', 'jensencompartmentexpts'),
    ('COMPARTMENTS Text-mining Protein Localization Evidence Scores', 'jensencompartmenttextmining'),
    ('CORUM Protein Complexes', 'corum'),
    ('COSMIC Cell Line Gene CNV Profiles', 'cosmiccnv'),
    ('COSMIC Cell Line Gene Mutation Profiles', 'cosmicmut'),
    ('CTD Gene-Chemical Interactions', 'ctdchemical'),
    ('CTD Gene-Disease Associations', 'ctddisease'),
    ('ClinVar SNP-Phenotype Associations', 'clinvar'),
    ('Combined Pathways Pathways', 'combinedpathways'),
    ('dbGAP Gene-Trait Associations', 'dbgap'),
    ('DEPOD Substrates of Phosphatases', 'depod'),
    ('DIP Protein-Protein Interactions', 'dip'),
    ('DISEASES Curated Gene-Disease Assocation Evidence Scores', 'jensendiseasecurated'),
    ('DISEASES Experimental Gene-Disease Assocation Evidence Scores', 'jensendiseaseexpts'),
    ('DISEASES Text-mining Gene-Disease Assocation Evidence Scores', 'jensendiseasetextmining'),
    ('DrugBank Drug Targets', 'drugbank'),
    ('ENCODE Histone Modification Site Profiles', 'encodehm'),
    ('ENCODE Transcription Factor Binding Site Profiles', 'encodetf'),
    ('ENCODE Transcription Factor Targets', 'encodetfppi'),
    ('ESCAPE Omics Signatures of Genes and Proteins for Stem Cells', 'escape'),
    ('GAD Gene-Disease Associations', 'gad'),
    ('GAD High Level Gene-Disease Associations', 'gadhighlevel'),
    ('GDSC Cell Line Gene Expression Profiles', 'gdsc'),
    ('GEO Signatures of Differentially Expressed Genes for Diseases', 'geodisease'),
    ('GEO Signatures of Differentially Expressed Genes for Gene Perturbations', 'geogene'),
    ('GEO Signatures of Differentially Expressed Genes for Kinase Perturbations', 'geokinase'),
    ('GEO Signatures of Differentially Expressed Genes for Small Molecules', 'geochemical'),
    ('GEO Signatures of Differentially Expressed Genes for Transcription Factor Perturbations', 'geotf'),
    ('GEO Signatures of Differentially Expressed Genes for Viral Infections', 'geovirus'),
    ('GO Biological Process Annotations', 'gobp'),
    ('GO Cellular Component Annotations', 'gocc'),
    ('GO Molecular Function Annotations', 'gomf'),
    ('GTEx Tissue Gene Expression Profiles', 'gtextissue'),
    ('GTEx Tissue Sample Gene Expression Profiles', 'gtexsample'),
    ('GTEx eQTL', 'gtexeqtl'),
    ('GWAS Catalog SNP-Phenotype Associations', 'gwascatalog'),
    ('GWASdb SNP-Disease Associations', 'gwasdbdisease'),
    ('GWASdb SNP-Phenotype Associations', 'gwasdbphenotype'),
    ('GeneRIF Biological Term Annotations', 'generif'),
    ('GeneSigDB Published Gene Signatures', 'genesigdb'),
    ('Guide to Pharmacology Chemical Ligands of Receptors', 'guidetopharmchemical'),
    ('Guide to Pharmacology Protein Ligands of Receptors', 'guidetopharmprotein'),
    ('HMDB Metabolites of Enzymes', 'hmdb'),
    ('HPA Cell Line Gene Expression Profiles', 'hpacelllines'),
    ('HPA Tissue Gene Expression Profiles', 'hpatissuesmrna'),
    ('HPA Tissue Protein Expression Profiles', 'hpatissuesprotein'),
    ('HPA Tissue Sample Gene Expression Profiles', 'hpasamples'),
    ('HPM Cell Type and Tissue Protein Expression Profiles', 'hpm'),
    ('HPO Gene-Disease Associations', 'hpo'),
    ('HPRD Protein-Protein Interactions', 'hprd'),
    ('Heiser et al., PNAS, 2011 Cell Line Gene Expression Profiles', 'heiser'),
    ('HuGE Navigator Gene-Phenotype Associations', 'hugenavigator'),
    ('Hub Proteins Protein-Protein Interactions', 'hubs'),
    ('HumanCyc Biomolecular Interactions', 'humancycppi'),
    ('HumanCyc Pathways', 'humancyc'),
    ('IntAct Biomolecular Interactions', 'intact'),
    ('InterPro Predicted Protein Domain Annotations', 'interpro'),
    ('JASPAR Predicted Transcription Factor Targets', 'jasparpwm'),
    ('KEA Substrates of Kinases', 'kea'),
    ('KEGG Biomolecular Interactions', 'keggppi'),
    ('KEGG Pathways', 'kegg'),
    ('Kinativ Kinase Inhibitor Bioactivity Profiles', 'kinativ'),
    ('KinomeScan Kinase Inhibitor Targets', 'kinomescan'),
    ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene CNV Profiles', 'klijncnv'),
    ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene Expression Profiles', 'klijnmrna'),
    ('Klijn et al., Nat. Biotechnol., 2015 Cell Line Gene Mutation Profiles', 'klijnmut'),
    ('LINCS L1000 CMAP Signatures of Differentially Expressed Genes for Small Molecules', 'lincscmapchemical'),
    ('LOCATE Curated Protein Localization Annotations', 'locate'),
    ('LOCATE Predicted Protein Localization Annotations', 'locatepredicted'),
    ('MPO Gene-Phenotype Associations', 'mgimpo'),
    ('MSigDB Cancer Gene Co-expression Modules', 'msigdbcomp'),
    ('MSigDB Signatures of Differentially Expressed Genes for Cancer Gene Perturbations', 'msigdbonc'),
    ('MiRTarBase microRNA Targets', 'mirtarbase'),
    ('MotifMap Predicted Transcription Factor Targets', 'motifmap'),
    ('NURSA Protein Complexes', 'nursa'),
    ('NURSA Protein-Protein Interactions', 'nursappi'),
    ('OMIM Gene-Disease Associations', 'omim'),
    ('PANTHER Biomolecular Interactions', 'pantherppi'),
    ('PANTHER Pathways', 'panther'),
    ('PID Biomolecular Interactions', 'pidppi'),
    ('PID Pathways', 'pid'),
    ('Pathway Commons Protein-Protein Interactions', 'pc'),
    ('PhosphoSitePlus Phosphosite-Disease Associations', 'phosphositeplusdisease'),
    ('PhosphoSitePlus Substrates of Kinases', 'phosphositeplus'),
    ('Phosphosite Textmining Biological Term Annotations', 'phosphositetextmining'),
    ('ProteomicsDB Cell Type and Tissue Protein Expression Profiles', 'proteomicsdb'),
    ('Reactome Biomolecular Interactions', 'reactomeppi'),
    ('Reactome Pathways', 'reactome'),
    ('Recon X Predicted Biomolecular Interactions', 'reconx'),
    ('Roadmap Epigenomics Cell and Tissue DNA Methylation Profiles', 'epigenomicsdnamethylation'),
    ('Roadmap Epigenomics Cell and Tissue Gene Expression Profiles', 'epigenomicsmrna'),
    ('Roadmap Epigenomics Histone Modification Site Profiles', 'epigenomicshm'),
    ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Drugs', 'silacdrug'),
    ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Gene Perturbations', 'silacgene'),
    ('SILAC Phosphoproteomics Signatures of Differentially Phosphorylated Proteins for Protein Ligands', 'silacligand'),
    ('TCGA Signatures of Differentially Expressed Genes for Tumors', 'tcga'),
    ('TISSUES Curated Tissue Protein Expression Evidence Scores', 'jensentissuecurated'),
    ('TISSUES Experimental Tissue Protein Expression Evidence Scores', 'jensentissueexpts'),
    ('TISSUES Text-mining Tissue Protein Expression Evidence Scores', 'jensentissuetextmining'),
    ('TRANSFAC Curated Transcription Factor Targets', 'transfac'),
    ('TRANSFAC Predicted Transcription Factor Targets', 'transfacpwm'),
    ('TargetScan Predicted Conserved microRNA Targets', 'targetscan'),
    ('TargetScan Predicted Nonconserved microRNA Targets', 'targetscannonconserved'),
    ('Virus MINT Protein-Viral Protein Interactions', 'virusmintppi'),
    ('Virus MINT Protein-Virus Interactions', 'virusmint'),
    ('Wikipathways Pathways', 'wikipathways'),
]

  selected_downloads = [
    'gene_attribute_matrix.txt.gz',
    'gene_attribute_edges.txt.gz',
    'gene_set_library_crisp.gmt.gz',
    'gene_set_library_up_crisp.gmt.gz',
    'gene_set_library_dn_crisp.gmt.gz',
    'attribute_set_library_crisp.gmt.gz',
    'attribute_set_library_up_crisp.gmt.gz',
    'attribute_set_library_dn_crisp.gmt.gz',
    'gene_similarity_matrix_cosine.txt.gz',
    'attribute_similarity_matrix_cosine.txt.gz',
    'gene_list_terms.txt.gz',
    'attribute_list_entries.txt.gz'
      ]

  for dataset, path in tqdm(selected_datasets):
      for downloadable in selected_downloads:
        url = 'https://maayanlab.cloud/static/hdfs/harmonizome/data/%s/%s' %\
          (path, downloadable)
        response = requests.get(url, stream=True)
          # Not every dataset has all downloadables.
        if response.status_code != 200:
          print("Download failed for %s %s, response status: %s" % (dataset,downloadable,response.status_code))
          continue   
        else:
            md5 = hashlib.md5()
            sha256 = hashlib.sha256()
            filesize = 0
            for chunk in response.iter_content(4096):
                sha256.update(chunk)
                md5.update(chunk)
                filesize += len(chunk)
        sha256_hash = sha256.hexdigest()
        md5_hash = md5.hexdigest()
          
        f = c2m2_level_1.file(
        id_namespace=ns.id,
        local_id=path+'_'+downloadable,
        project_id_namespace=ns.id,
        project_local_id='Harmonizome',
        filename=downloadable,
        size_in_bytes=filesize,
        persistent_id='https://maayanlab.cloud/static/hdfs/harmonizome/data/%s/%s' % (path,downloadable),
        sha256=sha256_hash,
        md5=md5_hash,
        )
        yield f


if __name__ == '__main__':
  if len(sys.argv) >= 2:
    outdir = sys.argv[1]
  else:
    logging.error("OUTDIR must be specified.")
    sys.exit()
  pkg = create_datapackage('C2M2_Level_1', convert_harmonizome_datasets_to_c2m2(), outdir)
  validate_datapackage(pkg)
  validate_id_namespace_name_uniqueness(pkg)