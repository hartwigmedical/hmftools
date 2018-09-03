# Knowledgebase-importer

The knowledgebase importer reads various knowledgebases and produces one consistent database holding:
 - Known Somatic SNV/Indel Hotspots
 - Known Somatic Fusions
 - Evidence for actionability (which variants may imply responsiveness or resistant to certain treatments?)
 
 The knowledgebase importer requires the full DOID definition (doid.owl) from  http://www.obofoundry.org/ontology/doid.html. The DOID definition is managed by http://www.disease-ontology.org.
 
 The following knowledgebases are integrated:
   - COSMIC for fusions (details to come...)
   - CGI for actionability and hotspots/fusions (https://www.cancergenomeinterpreter.org, and  https://doi.org/10.1101/140475)
   - OncoKB for actionability and hotspots/fusions (http://oncokb.org, and , and https://doi.org/10.1200/PO.17.00011)
   - CiViC for actionability and hotspots/fusions (https://civicdb.org, and https://www.ncbi.nlm.nih.gov/pubmed/28138153)
   - iClusion for (Dutch) clinical trials (https://iclusion.org)
   
 More specifically, the following files are required:
  - CGI biomarkers per variant (https://www.cancergenomeinterpreter.org/biomarkers): Download the folder and use file "cgi_biomarkers_per_variant.tsv".
  - CGI catalog of validated oncogenic mutations (https://www.cancergenomeinterpreter.org/mutations): Download the folder and use file "catalog_of_validated_oncogenic_mutations.tsv".
  - CiViC clinical evidence summaries (https://civicdb.org/releases): Use file "ClinicalEvidenceSummaries.tsv".
  - CiViC variant summaries (https://civicdb.org/releases): Use file "VariantSummaries.tsv".
  - OncoKB actionable variants (http://oncokb.org/#/dataAccess): Use link "actionable alterations".
  - OncoKB annotated variants (http://oncokb.org/#/dataAccess): Use link "all alterations".
  
  # iClusion
  
 iClusion provides personalized API access to their knowledgebase and has no publicly downloadable data available.
