SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS viccEntry;
DROP TABLE IF EXISTS gene;
DROP TABLE IF EXISTS geneIdentifier;
DROP TABLE IF EXISTS featureName;
DROP TABLE IF EXISTS tag;
DROP TABLE IF EXISTS devTag;
DROP TABLE IF EXISTS feature;
DROP TABLE IF EXISTS provenance;
DROP TABLE IF EXISTS synonym;
DROP TABLE IF EXISTS link;
DROP TABLE IF EXISTS sequenceOntology;
DROP TABLE IF EXISTS hierarchy;
DROP TABLE IF EXISTS association;
DROP TABLE IF EXISTS evidence;
DROP TABLE IF EXISTS evidenceInfo;
DROP TABLE IF EXISTS evidenceType;
DROP TABLE IF EXISTS publicationUrl;
DROP TABLE IF EXISTS phenotype;
DROP TABLE IF EXISTS phenotypeType;
DROP TABLE IF EXISTS environmentalContext;
DROP TABLE IF EXISTS approvedCountry;
DROP TABLE IF EXISTS taxonomy;
DROP TABLE IF EXISTS sage;
DROP TABLE IF EXISTS brca;

SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE viccEntry
(   id int NOT NULL AUTO_INCREMENT,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    geneName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE geneIdentifier
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    symbol varchar(255) NOT NULL,
    entrezId varchar(255) NOT NULL,
    ensemblGeneId varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);


CREATE TABLE featureName
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    nameOfFeature varchar(2500),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE tag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    tagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE devTag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    devTagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE feature
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    name varchar(1000),
    biomarkerType varchar(255),
    referenceName varchar(255),
    chromosome varchar(255),
    start varchar(255),
    end varchar(255),
    ref varchar(1000),
    alt varchar(255),
    provenanceRule varchar(255),
    geneSymbol varchar(255),
    entrezId varchar(255),
    description varchar(2000),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE provenance
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    provenanceName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE synonym
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    synonymName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE link
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    linkName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE sequenceOntology
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    soid varchar(255) NOT NULL,
    parentSoid varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    parentName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE hierarchy
(   id int NOT NULL AUTO_INCREMENT,
    sequenceOntologyId int NOT NULL,
    hierarchyName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (sequenceOntologyId) REFERENCES sequenceOntology(id)
);

CREATE TABLE association
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    variantName varchar(255),
    evidenceLevel varchar(255),
    evidenceLabel varchar(255),
    responseType varchar(255),
    drugLabel varchar(2000),
    sourceLink varchar(255),
    description varchar(2500) NOT NULL,
    oncogenic varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE evidence
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    description varchar(2000),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE evidenceInfo
(   id int NOT NULL AUTO_INCREMENT,
    evidenceId int NOT NULL,
    publication varchar(225) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id)
);

CREATE TABLE evidenceType
(   id int NOT NULL AUTO_INCREMENT,
    evidenceId int NOT NULL,
    sourceName varchar(225) NOT NULL,
    idEvidenceType varchar(225),
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id)
);

CREATE TABLE publicationUrl
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    urlOfPublication varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE phenotype
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    description varchar(255) NOT NULL,
    family varchar(255) NOT NULL,
    idPhenotype varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE phenotypeType
(   id int NOT NULL AUTO_INCREMENT,
    phenotypeId int NOT NULL,
    source varchar(255),
    term varchar(255) NOT NULL,
    idPhenotypeType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (phenotypeId) REFERENCES phenotype(id)
);

CREATE TABLE environmentalContext
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    term varchar(255),
    description varchar(255),
    source varchar(255),
    usanStem varchar(255),
    idEnvironmentalContext varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE approvedCountry
(   id int NOT NULL AUTO_INCREMENT,
    environmentContextId int NOT NULL,
    approvedCountryName varchar(225) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentContextId) REFERENCES environmentalContext(id)
);

CREATE TABLE taxonomy
(   id int NOT NULL AUTO_INCREMENT,
    environmentContextId int NOT NULL,
    kingdom varchar(225) NOT NULL,
    directParent varchar(225) NOT NULL,
    class varchar(225) NOT nULL,
    subClass varchar(225),
    superClass varchar(225) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentContextId) REFERENCES environmentalContext(id)
);

CREATE TABLE sage
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    entrezId varchar(255) NOT NULL,
    clinicalManifestation varchar(255) NOT NULL,
    publicationUrl varchar(255) NOT NULL,
    germlineOrSomatic varchar(255) NOT NULL,
    evidenceLabel varchar(255) NOT NULL,
    drugLabel varchar(255) NOT NULL,
    responseType varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

/*
CREATE TABLE brca
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    Variant_frequency_LOVD varchar(255) NOT NULL,
    Allele_frequency_FIN_ExAC varchar(255) NOT NULL,
    ClinVarAccession_ENIGMA varchar(255) NOT NULL,
    Homozygous_count_AFR_ExAC varchar(255) NOT NULL,
    BX_ID_ExAC varchar(255) NOT NULL,
    Variant_in_LOVD varchar(255) NOT NULL,
    Allele_frequency_AFR_ExAC varchar(255) NOT NULL,
    DBID_LOVD varchar(255) NOT NULL,
    Chr varchar(255) NOT NULL,
    BX_ID_ENIGMA varchar(255) NOT NULL,
    Co_occurrence_LR_exLOVD varchar(255) NOT NULL,
    Homozygous_count_EAS_ExAC varchar(255) NOT NULL,
    Submitter_ClinVar varchar(1500) NOT NULL,
    Allele_frequency_EAS_ExAC varchar(255) NOT NULL,
    Hg37_End varchar(255) NOT NULL,
    Submitters_LOVD varchar(1000) NOT NULL,
    Clinical_classification_BIC varchar(255) NOT NULL,
    Homozygous_count_NFE_ExAC varchar(255) NOT NULL,
    Allele_count_SAS_ExAC varchar(255) NOT NULL,
    Method_ClinVar varchar(255) NOT NULL,
    Allele_count_NFE_ExAC varchar(255) NOT NULL,
    Pathogenicity_all varchar(255) NOT NULL,
    Germline_or_Somatic_BIC varchar(255) NOT NULL,
    Homozygous_count_SAS_ExAC varchar(255) NOT NULL,
    BIC_Nomenclature varchar(255) NOT NULL,
    Assertion_method_ENIGMA varchar(255) NOT NULL,
    Literature_source_exLOVD varchar(255) NOT NULL,
    Change_Type_id varchar(255) NOT NULL,
    Collection_method_ENIGMA varchar(255) NOT NULL,
    Sum_family_LR_exLOVD varchar(255) NOT NULL,
    HGVS_cDNA_LOVD varchar(255) NOT NULL,
    Homozygous_count_FIN_ExAC varchar(255) NOT NULL,
    EAS_Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    Ethnicity_BIC varchar(2000) NOT NULL,
    Individuals_LOVD varchar(255) NOT NULL,
    Variant_in_ExAC varchar(255) NOT NULL,
    URL_ENIGMA varchar(255) NOT NULL,
    Allele_Origin_ClinVar varchar(255) NOT NULL,
    Allele_frequency_AMR_ExAC varchar(255) NOT NULL,
    Variant_in_1000_Genomes varchar(255) NOT NULL,
    AFR_Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    BX_ID_exLOVD varchar(255) NOT NULL,
    Source varchar(255) NOT NULL,
    Condition_ID_value_ENIGMA varchar(255) NOT NULL,
    HGVS_Protein varchar(255) NOT NULL,
    Ref varchar(255) NOT NULL,
    Allele_number_AFR_ExAC varchar(255) NOT NULL,
    Allele_count_AFR_ExAC varchar(255) NOT NULL,
    BX_ID_LOVD varchar(255) NOT NULL,
    Synonyms varchar(1000) NOT NULL,
    Gene_Symbol varchar(255) NOT NULL,
    Comment_on_clinical_significance_ENIGMA varchar(500) NOT NULL,
    Missense_analysis_prior_probability_exLOVD varchar(255) NOT NULL,
    Allele_number_FIN_ExAC varchar(255) NOT NULL,
    Posterior_probability_exLOVD varchar(255) NOT NULL,
    Polyphen_Score varchar(255) NOT NULL,
    Reference_Sequence varchar(255) NOT NULL,
    Allele_count_EAS_ExAC varchar(255) NOT NULL,
    Hg38_End varchar(255) NOT NULL,
    HGVS_cDNA varchar(255) NOT NULL,
    Functional_analysis_technique_LOVD varchar(255) NOT NULL,
    SAS_Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    RNA_LOVD varchar(255) NOT NULL,
    Combined_prior_probablility_exLOVD varchar(255) NOT NULL,
    BX_ID_ClinVar varchar(255) NOT NULL,
    IARC_class_exLOVD varchar(255) NOT NULL,
    BX_ID_BIC varchar(12500) NOT NULL,
    Sift_Prediction varchar(255) NOT NULL,
    Allele_number_NFE_ExAC varchar(255) NOT NULL,
    Allele_origin_ENIGMA varchar(255) NOT NULL,
    Allele_number_OTH_ExAC varchar(255) NOT NULL,
    Hg36_End varchar(255) NOT NULL,
    Allele_frequency_SAS_ExAC varchar(255) NOT NULL,
    Date_Last_Updated_ClinVar varchar(500) NOT NULL,
    Allele_number_EAS_ExAC varchar(255) NOT NULL,
    Allele_frequency_OTH_ExAC varchar(255) NOT NULL,
    Source_URL varchar(2000) NOT NULL,
    SCV_ClinVar varchar(500) NOT NULL,
    Pathogenicity_expert varchar(255) NOT NULL,
    Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    Functional_analysis_result_LOVD varchar(255) NOT NULL,
    AMR_Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    Variant_in_ESP varchar(255) NOT NULL,
    Variant_in_BIC varchar(255) NOT NULL,
    Clinical_significance_ENIGMA varchar(255) NOT NULL,
    Max_Allele_Frequency varchar(255) NOT NULL,
    Allele_count_AMR_ExAC varchar(255) NOT NULL,
    Variant_in_ENIGMA varchar(255) NOT NULL,
    BX_ID_ESP varchar(255) NOT NULL,
    Patient_nationality_BIC varchar(500) NOT NULL,
    BX_ID_1000_Genomes varchar(255) NOT NULL,
    Genomic_Coordinate_hg37 varchar(255) NOT NULL,
    Genomic_Coordinate_hg36 varchar(255) NOT NULL,
    EUR_Allele_frequency_1000_Genomes varchar(255) NOT NULL,
    Number_of_family_member_carrying_mutation_BIC varchar(255) NOT NULL,
    Segregation_LR_exLOVD varchar(255) NOT NULL,
    Allele_Frequency varchar(255) NOT NULL,
    Minor_allele_frequency_percent_ESP varchar(255) NOT NULL,
    Allele_frequency_ExAC varchar(255) NOT NULL,
    Mutation_type_BIC varchar(255) NOT NULL,
    Assertion_method_citation_ENIGMA varchar(255) NOT NULL,
    Condition_ID_type_ENIGMA varchar(255) NOT NULL,
    Allele_count_OTH_ExAC varchar(255) NOT NULL,
    HGVS_protein_LOVD varchar(255) NOT NULL,
    Variant_in_ClinVar varchar(255) NOT NULL,
    Clinical_importance_BIC varchar(255) NOT NULL,
    Discordant varchar(255) NOT NULL,
    Allele_count_FIN_ExAC varchar(255) NOT NULL,
    Condition_category_ENIGMA varchar(255) NOT NULL,
    Allele_Frequency_ESP varchar(255) NOT NULL,
    Homozygous_count_OTH_ExAC varchar(255) NOT NULL,
    Genetic_origin_LOVD varchar(255) NOT NULL,
    Homozygous_count_AMR_ExAC varchar(255) NOT NULL,
    Clinical_Significance_ClinVar varchar(255) NOT NULL,
    AA_Allele_Frequency_ESP varchar(255) NOT NULL,
    Protein_Change varchar(255) NOT NULL,
    Variant_in_exLOVD varchar(255) NOT NULL,
    EA_Allele_Frequency_ESP varchar(255) NOT NULL,
    HGVS_RNA varchar(255) NOT NULL,
    Clinical_significance_citations_ENIGMA varchar(255) NOT NULL,
    Variant_effect_LOVD varchar(255) NOT NULL,
    Polyphen_Prediction varchar(255) NOT NULL,
    Data_Release_id varchar(255) NOT NULL,
    Hg37_Start varchar(255) NOT NULL,
    Hg36_Start varchar(255) NOT NULL,
    Sift_Score varchar(255) NOT NULL,
    Genomic_Coordinate_hg38 varchar(255) NOT NULL,
    Alt varchar(255) NOT NULL,
    Literature_citation_BIC varchar(1000) NOT NULL,
    Variant_haplotype_LOVD varchar(255) NOT NULL,
    Allele_frequency_NFE_ExAC varchar(255) NOT NULL,
    Hg38_Start varchar(255) NOT NULL,
    Pos varchar(255) NOT NULL,
    Date_last_evaluated_ENIGMA varchar(255) NOT NULL,
    Allele_number_SAS_ExAC varchar(255) NOT NULL,
    Allele_number_AMR_ExAC varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);*/
