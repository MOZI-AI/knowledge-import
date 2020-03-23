
Bio-Atomspace import scripts
============================

The Atomspace of the gene annotation service (MOZI-AI) integrates the follwoing 
list of biology datasets (of Homo sapiens species).

1. The Biogrid: BioGRID is a repository for Interaction Datasets. 
   - biogrid.py imports Biogrid interaction of Human genes
   - gene2proteinMapping.py imports the uniprot proteins and the genes expressing them
   - uniprotIDmap.R generates the gene-protein mapping file

2. Reactome pathway: The complete list of pathways and hierarchial relationship among them.
   - Use reactome.py to get that in atomese version

   The three Physical Entity (PE) Identifier mapping files 
   Physical Entity (PE) Identifier mapping.py imports the following
	- Mapping Uniprots to Pathways
	- Mapping ChEbi to Pathways and
	- NCBI Genes to pathways

3. Small molecule Pathway database (SMPDB)

   The Metabolite names linked to SMPDB pathways and Protein names linked to SMPDB pathways

4. Gene ontology database

   The Genes and their ontology GO (classes used to describe gene function
   and relationships betweeen these classes)

5. STRING Protein-Protein Interaction Networks Functional Enrichment Analysis

The imported Atomese version of the datasets can be found https://mozi.ai/datasets/

NOTE: For expermenting only on a gene-level, the following scripts generate a reduced version of the data

- biogrid.py
- GO_Annotation_scm.py
- SMPDB_pathway.py
