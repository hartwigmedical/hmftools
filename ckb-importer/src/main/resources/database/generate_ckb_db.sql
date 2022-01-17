SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS ckbEntry;
CREATE TABLE ckbEntry
(   id int NOT NULL AUTO_INCREMENT,
    ckbProfileId int NOT NULL UNIQUE,
    createDate DATE NOT NULL,
    updateDate DATE NOT NULL,
    profileName varchar(250) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS variant;
CREATE TABLE variant
(   id int NOT NULL AUTO_INCREMENT,
    ckbEntryId int NOT NULL,
    ckbVariantId int NOT NULL,
    createDate DATE NOT NULL,
    updateDate DATE NOT NULL,
    fullName varchar(50) NOT NULL,
    variant varchar(50) NOT NULL,
    impact varchar(50),
    proteinEffect varchar(50),
    type varchar(50),
    description varchar(2500),
    PRIMARY KEY (id),
    FOREIGN KEY (ckbEntryId) REFERENCES ckbEntry(id)
);

DROP TABLE IF EXISTS variantReference;
CREATE TABLE variantReference
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    ckbReferenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
);

DROP TABLE IF EXISTS gene;
CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    ckbGeneId int NOT NULL,
    createDate DATE NOT NULL,
    updateDate DATE,
    geneSymbol varchar(50) NOT NULL,
    geneRole varchar(250) NOT NULL,
    entrezId varchar(50),
    chromosome varchar(50),
    mapLocation varchar(50),
    canonicalTranscript varchar(50),
    description varchar(2500),
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

DROP TABLE IF EXISTS geneReference;
CREATE TABLE geneReference
(   id int NOT NULL AUTO_INCREMENT,
    geneId int NOT NULL,
    ckbReferenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (geneId) REFERENCES gene(id)
);

DROP TABLE IF EXISTS transcriptCoordinate;
CREATE TABLE transcriptCoordinate
(   id int NOT NULL AUTO_INCREMENT,
    variantId int NOT NULL,
    isReferenceTranscriptCoordinate BOOLEAN NOT NULL,
    transcript varchar(50) NOT NULL,
    gDna varchar(500) NOT NULL,
    cDna varchar(500) NOT NULL,
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
    variantPath varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (variantId) REFERENCES variant(id)
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

DROP TABLE IF EXISTS evidence;
CREATE TABLE evidence
(   id int NOT NULL AUTO_INCREMENT,
    ckbEntryId int NOT NULL,
    ckbEvidenceId int NOT NULL,
    responseType varchar(50) NOT NULL,
    evidenceType varchar(50) NOT NULL,
    efficacyEvidence varchar(5000) NOT NULL,
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
    ckbReferenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
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
    createDate DATE NOT NULL,
    updateDate DATE,
    therapyName varchar(500) NOT NULL,
    description varchar(2500),
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
    synonym varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS therapyReference;
CREATE TABLE therapyReference
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    ckbReferenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS globalApprovalStatus;
CREATE TABLE globalApprovalStatus
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    ckbGlobalApprovalStatusId int NOT NULL,
    ckbProfileId int NOT NULL,
    ckbIndicationId int NOT NULL,
    approvalStatus varchar(50) NOT NULL,
    approvalAuthority varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS drug;
CREATE TABLE drug
(   id int NOT NULL AUTO_INCREMENT,
    therapyId int NOT NULL,
    ckbDrugId int NOT NULL,
    createDate DATE NOT NULL,
    drugName varchar(50) NOT NULL,
    tradeName varchar(50),
    casRegistryNum varchar(50),
    ncitId varchar(50),
    description varchar(2500),
    PRIMARY KEY (id),
    FOREIGN KEY (therapyId) REFERENCES therapy(id)
);

DROP TABLE IF EXISTS drugClass;
CREATE TABLE drugClass
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    ckbDrugClassId int NOT NULL,
    createDate DATE NOT NULL,
    drugClass varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugTerm;
CREATE TABLE drugTerm
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    term varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugSynonym;
CREATE TABLE drugSynonym
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    synonym varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS drugReference;
CREATE TABLE drugReference
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    ckbReferenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

DROP TABLE IF EXISTS indication;
CREATE TABLE indication
(   id int NOT NULL AUTO_INCREMENT,
    ckbIndicationId int NOT NULL,
    name varchar(250) NOT NULL,
    source varchar(50) NOT NULL,
    definition varchar(5000),
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
    title varchar(500) NOT NULL,
    phase varchar(50) NOT NULL,
    recruitment varchar(50) NOT NULL,
    gender varchar(50),
    sponsors varchar(250),
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

DROP TABLE IF EXISTS ageGroup;
CREATE TABLE ageGroup
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ageGroup varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS variantRequirementDetail;
CREATE TABLE variantRequirementDetail
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ckbProfileId int NOT NULL,
    requirementType varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS location;
CREATE TABLE location
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    nctId varchar(50) NOT NULL,
    status varchar(50),
    facility varchar(500),
    city varchar(50) NOT NULL,
    state varchar(50),
    zip varchar(50),
    country varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS contact;
CREATE TABLE contact
(   id int NOT NULL AUTO_INCREMENT,
    locationId int NOT NULL,
    name varchar(250),
    email varchar(250),
    phone varchar(50),
    phoneExt varchar(50),
    role varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (locationId) REFERENCES location(id)
);

DROP TABLE IF EXISTS treatmentApproach;
CREATE TABLE treatmentApproach
(   id int NOT NULL AUTO_INCREMENT,
    treatmentApproachId int NOT NULL,
    createDate DATE NOT NULL,
    updateDate DATE,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS treatmentApproachEvidence;
CREATE TABLE treatmentApproachEvidence
(   evidenceId int NOT NULL,
    treatmentApproachEvidenceId int NOT NULL,
    PRIMARY KEY (evidenceId, treatmentApproachEvidenceId),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id),
    FOREIGN KEY (treatmentApproachEvidenceId) REFERENCES treatmentApproach(id)
);

DROP TABLE IF EXISTS treatmentApproachDrugClass;
CREATE TABLE treatmentApproachDrugClass
(   id int NOT NULL AUTO_INCREMENT,
    treatmentApproachId int NOT NULL,
    drugClassId int NOT NULL,
    createDate DATE NOT NULL,
    drugClass varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentApproachId) REFERENCES treatmentApproach(id)
);

DROP TABLE IF EXISTS treatmentApproachReference;
CREATE TABLE treatmentApproachReference
(   id int NOT NULL AUTO_INCREMENT,
    treatmentApproachId int NOT NULL,
    referenceId int NOT NULL,
    pubMedId varchar(50),
    title varchar(500),
    abstractText TEXT,
    url varchar(250),
    journal varchar(500),
    authors varchar(5000),
    volume varchar(50),
    issue varchar(50),
    date varchar(50),
    year varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentApproachId) REFERENCES treatmentApproach(id)
);

SET FOREIGN_KEY_CHECKS = 1;