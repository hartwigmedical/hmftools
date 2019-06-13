SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS viccEntry;
CREATE TABLE viccEntry
(   id int NOT NULL AUTO_INCREMENT,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS gene;
CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    geneName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS geneIdentifier;
CREATE TABLE geneIdentifier
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    symbol varchar(255) NOT NULL,
    entrezId varchar(255) NOT NULL,
    ensemblGeneId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS featureName;
CREATE TABLE featureName
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    featureName varchar(2500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS tag;
CREATE TABLE tag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    tagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS devTag;
CREATE TABLE devTag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    devTagName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS feature;
CREATE TABLE feature
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    name varchar(255),
    biomarkerType varchar(255),
    referenceName varchar(255),
    chromosome varchar(255),
    start varchar(255),
    end varchar(255),
    ref varchar(255),
    alt varchar(255),
    provenanceRule varchar(255),
    geneSymbol varchar(255),
    entrezId varchar(255),
    description varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS provenance;
CREATE TABLE provenance
(   id int NOT NULL AUTO_INCREMENT,
    featureEntryId int NOT NULL,
    provenance varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureEntryId) REFERENCES feature(id)
);

DROP TABLE IF EXISTS synonyms;
CREATE TABLE synonyms
(   id int NOT NULL AUTO_INCREMENT,
    featureEntryId int NOT NULL,
    synonyms varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureEntryId) REFERENCES feature(id)
);

DROP TABLE IF EXISTS links;
CREATE TABLE links
(   id int NOT NULL AUTO_INCREMENT,
    featureEntryId int NOT NULL,
    links varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureEntryId) REFERENCES feature(id)
);

DROP TABLE IF EXISTS sequenceOntology;
CREATE TABLE sequenceOntology
(   id int NOT NULL AUTO_INCREMENT,
    featureEntryId int NOT NULL,
    soid varchar(255),
    parentSoid varchar(255),
    name varchar(255),
    parentName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (featureEntryId) REFERENCES feature(id)
);

DROP TABLE IF EXISTS hierarchy;
CREATE TABLE hierarchy
(   id int NOT NULL AUTO_INCREMENT,
    sequenceOntologyEntryId int NOT NULL,
    hierarchy varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (sequenceOntologyEntryId) REFERENCES sequenceOntology(id)
);

DROP TABLE IF EXISTS association;
CREATE TABLE association
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    variantName varchar(255),
    evidenceLevel varchar(255),
    evidenceLabel varchar(255),
    responseType varchar(255),
    drugLabels varchar(255),
    sourceLink varchar(255),
    description varchar(500),
    oncogenic varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS evidence;
CREATE TABLE evidence
(   id int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    description varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationEntryId) REFERENCES association(id)
);

DROP TABLE IF EXISTS evidenceInfo;
CREATE TABLE evidenceInfo
(   id int NOT NULL AUTO_INCREMENT,
    evidenceEntryId int NOT NULL,
    publications varchar(225),
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceEntryId) REFERENCES evidence(id)
);

DROP TABLE IF EXISTS evidenceType;
CREATE TABLE evidenceType
(   id int NOT NULL AUTO_INCREMENT,
    evidenceEntryId int NOT NULL,
    sourceName varchar(225),
    idEvidenceType varchar(225),
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceEntryId) REFERENCES evidence(id)
);

DROP TABLE IF EXISTS publicationUrl;
CREATE TABLE publicationUrl
(   id int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    publicationUrls varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationEntryId) REFERENCES association(id)
);

DROP TABLE IF EXISTS phenotype;
CREATE TABLE phenotype
(   id int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    description varchar(255),
    family varchar(255),
    idPhenotype varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationEntryId) REFERENCES association(id)
);

DROP TABLE IF EXISTS phenotypeType;
CREATE TABLE phenotypeType
(   id int NOT NULL AUTO_INCREMENT,
    phenotypeEntryId int NOT NULL,
    source varchar(255),
    term varchar(255),
    idPhenotypeType varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (phenotypeEntryId) REFERENCES phenotype(id)
);

DROP TABLE IF EXISTS environmentalContext;
CREATE TABLE environmentalContext
(   id int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    term varchar(255),
    description varchar(255),
    source varchar(255),
    usanStem varchar(255),
    idEnvironmentalContexts varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (associationEntryId) REFERENCES association(id)
);

DROP TABLE IF EXISTS approvedCountries;
CREATE TABLE approvedCountries
(   id int NOT NULL AUTO_INCREMENT,
    environmentalContextsEntryId int NOT NULL,
    approvedCountries varchar(225),
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextsEntryId) REFERENCES environmentalContext(id)
);

DROP TABLE IF EXISTS taxonomy;
CREATE TABLE taxonomy
(   id int NOT NULL AUTO_INCREMENT,
    environmentalContextsEntryId int NOT NULL,
    kingdom varchar(225),
    directParent varchar(225),
    class varchar(225),
    subClass varchar(225),
    superClass varchar(225),
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextsEntryId) REFERENCES environmentalContext(id)
);

SET FOREIGN_KEY_CHECKS = 1;