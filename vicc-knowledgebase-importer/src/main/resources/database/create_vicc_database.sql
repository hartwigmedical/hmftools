SET FOREIGN_KEY_CHECKS = 0;

-- TODO: Remove per 1st of july.
DROP TABLE IF EXISTS features;

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
    environmentalContextsId int NOT NULL,
    approvedCountryName varchar(225) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextsId) REFERENCES environmentalContext(id)
);

CREATE TABLE taxonomy
(   id int NOT NULL AUTO_INCREMENT,
    environmentalContextsId int NOT NULL,
    kingdom varchar(225) NOT NULL,
    directParent varchar(225) NOT NULL,
    class varchar(225) NOT nULL,
    subClass varchar(225),
    superClass varchar(225) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextsId) REFERENCES environmentalContext(id)
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