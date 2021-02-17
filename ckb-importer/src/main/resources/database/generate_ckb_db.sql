SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS ckbEntry;
CREATE TABLE ckbEntry
(   id int NOT NULL AUTO_INCREMENT,
    ckbProfileId int NOT NULL UNIQUE,
    profileName varchar(50) NOT NULL,
    createDate DATE,
    updateDate DATE,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS variant;
CREATE TABLE variant
(   id int NOT NULL AUTO_INCREMENT,
    ckbEntryId int NOT NULL,
    ckbVariantId int NOT NULL,
    fullName varchar(50) NOT NULL,
    impact varchar(50),
    proteinEffect varchar(50),
    type varchar(50),
    variant varchar(50) NOT NULL,
    createDate DATE,
    updateDate DATE,
    PRIMARY KEY (id),
    FOREIGN KEY (ckbEntryId) REFERENCES ckbEntry(id)
);

DROP TABLE IF EXISTS gene;
CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    ckbGeneId int NOT NULL,
    geneSymbol varchar(50) NOT NULL,
    entrezId varchar(50),
    chromosome varchar(50),
    mapLocation varchar(50),
    canonicalTranscript varchar(50),
    geneRole varchar(50) NOT NULL,
    createDate DATE,
    updateDate DATE,
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS geneTerm;
CREATE TABLE geneTerm
(   id int NOT NULL AUTO_INCREMENT,
    geneId int NOT NULL,
    term varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (geneId) REFERENCES gene(id)
);

DROP TABLE IF EXISTS geneSynonym;
CREATE TABLE geneSynonym
(   id int NOT NULL AUTO_INCREMENT,
    geneId int NOT NULL,
    synonym varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (geneId) REFERENCES gene(id)
);

DROP TABLE IF EXISTS geneDescription;
CREATE TABLE geneDescription
(   id int NOT NULL AUTO_INCREMENT,
    geneId int NOT NULL,
    description varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (geneId) REFERENCES gene(id)
);

DROP TABLE IF EXISTS geneDescriptionReference;
CREATE TABLE geneDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    geneDescriptionId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (geneDescriptionId) REFERENCES geneDescription(id)
);

DROP TABLE IF EXISTS variantDescription;
CREATE TABLE variantDescription
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    description varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS variantDescriptionReference;
CREATE TABLE variantDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    variantDescriptionId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (variantDescriptionId) REFERENCES variantDescription(id)
);

DROP TABLE IF EXISTS referenceTranscriptCoordinate;
CREATE TABLE referenceTranscriptCoordinate
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    isReferenceTranscriptCoordinate BOOLEAN NOT NULL,
    transcript varchar(50) NOT NULL,
    gDna varchar(50) NOT NULL,
    cDna varchar(50) NOT NULL,
    protein varchar(50) NOT NULL,
    sourceDb varchar(50) NOT NULL,
    refGenomeBuild varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS categoryVariantPath;
CREATE TABLE categoryVariantPath
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    variantPath varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS categoryVariant;
CREATE TABLE categoryVariant
(   id int NOT NULL AUTO_INCREMENT,
    categoryVariantPathId int NOT NULL,
    ckbVariantId int NOT NULL,
    fullName varchar(50) NOT NULL,
    impact varchar(50),
    proteinEffect varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (categoryVariantPathId) REFERENCES categoryVariantPath(id)
);

DROP TABLE IF EXISTS memberVariant;
CREATE TABLE memberVariant
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    ckbVariantId int NOT NULL,
    fullName varchar(50) NOT NULL,
    impact varchar(50),
    proteinEffect varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS memberVariantDescription;
CREATE TABLE memberVariantDescription
(   id int NOT NULL AUTO_INCREMENT,
    memberVariantId int NOT NULL,
    description varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (memberVariantId) REFERENCES memberVariant(id)
);

DROP TABLE IF EXISTS memberVariantDescriptionReference;
CREATE TABLE memberVariantDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    memberVariantDescriptionId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (memberVariantDescriptionId) REFERENCES memberVariantDescription(id)
);

DROP TABLE IF EXISTS evidence;
CREATE TABLE evidence
(   id int NOT NULL AUTO_INCREMENT,
    ckbEntryId int NOT NULL,
    ckbEvidenceId int NOT NULL,
    responseType varchar(50) NOT NULL,
    evidenceType varchar(50) NOT NULL,
    efficacyEvidence varchar(50) NOT NULL,
    approvalStatus varchar(50) NOT NULL,
    ampCapAscoEvidenceLevel varchar(50) NOT NULL,
    ampCapAscoInferredTier varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (ckbEntryId) REFERENCES ckbEntry(id)
);

DROP TABLE IF EXISTS evidenceReference;
CREATE TABLE evidenceReference
(   id int NOT NULL AUTO_INCREMENT,
    evidenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id)
);

DROP TABLE IF EXISTS therapy;
CREATE TABLE therapy
(   id int NOT NULL AUTO_INCREMENT,
    ckbTherapyId int NOT NULL,
    createDate DATE,
    updateDate DATE,
    therapyName varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS therapyEvidence;
CREATE TABLE therapyEvidence
(   evidenceId int NOT NULL,
    therapyId int NOT NULL,
    PRIMARY KEY (evidenceId, therapyId),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS therapySynonym;
CREATE TABLE therapySynonym
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    synonym varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS therapyDescription;
CREATE TABLE therapyDescription
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    description varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS therapyDescriptionReference;
CREATE TABLE therapyDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    therapyDescriptionId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (therapyDescriptionId) REFERENCES therapyDescription(id)
);

DROP TABLE IF EXISTS globalTherapyApprovalStatus;
CREATE TABLE globalTherapyApprovalStatus
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    ckbGlobalTherapyApprovalStatusId int NOT NULL,
    ckbProfileId int NOT NULL,
    ckbTherapyId int NOT NULL,
    ckbIndicationId varchar(50) NOT NULL,
    approvalStatus varchar(50) NOT NULL,
    approvalAuthority varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS drug;
CREATE TABLE drug
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    ckbDrugId int NOT NULL,
    createDate DATE,
    drugName varchar(50) NOT NULL,
    tradeName varchar(50),
    casRegistryNum varchar(50),
    ncitId varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS drugClass;
CREATE TABLE drugClass
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    ckbDrugClassId int NOT NULL,
    createDate DATE,
    drugClass varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugTerm;
CREATE TABLE drugTerm
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    term varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugSynonym;
CREATE TABLE drugSynonym
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    synonym varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugDescription;
CREATE TABLE drugDescription
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    description varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugDescriptionReference;
CREATE TABLE drugDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    drugDescriptionId int NOT NULL,
    pubMedId varchar(50),
    title varchar(50),
    abstractText varchar(50),
    url varchar(50),
    authors varchar(50),
    journal varchar(50),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (drugDescriptionId) REFERENCES drugDescription(id)
);

DROP TABLE IF EXISTS indication;
CREATE TABLE indication
(   id int NOT NULL AUTO_INCREMENT,
    ckbIndicationId int NOT NULL,
    name varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    definition varchar(50),
    currentPreferredTerm varchar(50),
    lastUpdateDateFromDO DATE,
    termId varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS indicationEvidence;
CREATE TABLE indicationEvidence
(   evidenceId int NOT NULL,
    indicationId int NOT NULL,
    PRIMARY KEY (evidenceId, indicationId),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id),
    FOREIGN KEY (indicationId) REFERENCES indication(id)
);

DROP TABLE IF EXISTS indicationAltId;
CREATE TABLE indicationAltId
(   id int NOT NULL AUTO_INCREMENT,
    indicationId int NOT NULL,
    altId varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (indicationId) REFERENCES indication(id)
);

DROP TABLE IF EXISTS clinicalTrial;
CREATE TABLE clinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    ckbEntryId int NOT NULL,
    updateDate DATE,
    nctId varchar(50) NOT NULL,
    title varchar(50) NOT NULL,
    phase varchar(50) NOT NULL,
    recruitment varchar(50) NOT NULL,
    gender varchar(50),
    sponsor varchar(50),
    variantRequirement varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (ckbEntryId) REFERENCES ckbEntry(id)
);

DROP TABLE IF EXISTS therapyClinicalTrial;
CREATE TABLE therapyClinicalTrial
(   clinicalTrialId int NOT NULL,
    therapyId int NOT NULL,
    PRIMARY KEY (clinicalTrialId, therapyId),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS indicationClinicalTrial;
CREATE TABLE indicationClinicalTrial
(   clinicalTrialId int NOT NULL,
    indicationId int NOT NULL,
    PRIMARY KEY (clinicalTrialId, indicationId),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id),
    FOREIGN KEY (indicationId) REFERENCES indication(id)
);

DROP TABLE IF EXISTS clinicalTrialAgeGroup;
CREATE TABLE clinicalTrialAgeGroup
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ageGroup varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetail;
CREATE TABLE clinicalTrialVariantRequirementDetail
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ckbProfileId int NOT NULL,
    requirementType varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialLocation;
CREATE TABLE clinicalTrialLocation
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    nctId varchar(50) NOT NULL,
    facility varchar(50),
    city varchar(50) NOT NULL,
    country varchar(50) NOT NULL,
    status varchar(50),
    state varchar(50),
    zip varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialContact;
CREATE TABLE clinicalTrialContact
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialLocationId int NOT NULL,
    name varchar(50),
    email varchar(50),
    phone varchar(50),
    phoneExt varchar(50),
    role varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialLocationId) REFERENCES clinicalTrialLocation(id)
);

SET FOREIGN_KEY_CHECKS = 1;