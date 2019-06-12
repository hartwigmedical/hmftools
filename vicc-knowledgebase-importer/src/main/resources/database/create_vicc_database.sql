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
(   idFeature int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    name varchar(255),
    biomarkerType varchar(255),
    referenceName varchar(255),
    chromosome varchar(255),
    start varchar(255),
    end varchar(255),
    ref varchar(255),
    alt varchar(255),
    provenance varchar(255),
    provenanceRule varchar(255),
    geneSymbol varchar(255),
    synonyms varchar(255),
    entrezId varchar(255),
    links varchar(255),
    description varchar(255),
    PRIMARY KEY (idFeature),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS sequenceOntology;
CREATE TABLE sequenceOntology
(   id int NOT NULL AUTO_INCREMENT,
    featureEntryId int NOT NULL,
    hierarchy varchar(255),
    soid varchar(255) NOT NULL,
    parentSoid varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    parentName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureEntryId) REFERENCES feature(idFeature)
);

DROP TABLE IF EXISTS association;
CREATE TABLE association
(   idAssociation int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    variantName varchar(255),
    evidenceLevel varchar(255),
    evidenceLabel varchar(255),
    responseType varchar(255),
    drugLabels varchar(255),
    sourceLink varchar(255),
    publicationUrls varchar(255),
    description varchar(255) NOT NULL,
    oncogenic varchar(255),
    PRIMARY KEY (idAssociation),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

DROP TABLE IF EXISTS evidence;
CREATE TABLE evidence
(   idEvidence int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    evidenceInfoPublications varchar(255),
    evidenceTypeSourceName varchar(255),
    evidenceTypeId varchar(255),
    description varchar(255),
    PRIMARY KEY (idEvidence),
    FOREIGN KEY (associationEntryId) REFERENCES association(idAssociation)
);

DROP TABLE IF EXISTS phenotype;
CREATE TABLE phenotype
(   idPhenotype int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    phenotypeTypeSource varchar(255),
    phenotypeTypeterm varchar(255),
    phenotypeTypeid varchar(255),
    description varchar(255),
    family varchar(255),
    id varchar(255),
    PRIMARY KEY (idPhenotype),
    FOREIGN KEY (associationEntryId) REFERENCES association(idAssociation)
);

DROP TABLE IF EXISTS environmentalContext;
CREATE TABLE environmentalContext
(   idEnvironmentalContext int NOT NULL AUTO_INCREMENT,
    associationEntryId int NOT NULL,
    term varchar(255),
    description varchar(255),
    taxonomyKingdom varchar(255),
    taxonomyDirectParent varchar(255),
    taxonomyClass varchar(255),
    taxonomySubClass varchar(255),
    taxonomySuperClass varchar(255),
    source varchar(255),
    usanStem varchar(255),
    approvedCountries varchar(255),
    id varchar(255),
    PRIMARY KEY (idEnvironmentalContext),
    FOREIGN KEY (associationEntryId) REFERENCES association(idAssociation)
);


SET FOREIGN_KEY_CHECKS = 1;