SET FOREIGN_KEY_CHECKS = 0;

-- TODO Remove per 1st of March 2021
DROP TABLE IF EXISTS ckbEntry;

DROP TABLE IF EXISTS clinicalTrialContact;
DROP TABLE IF EXISTS clinicalTrialLocation;
DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetailMolecularProfile;
DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetail;
DROP TABLE IF EXISTS clinicalTrialIndication;
DROP TABLE IF EXISTS clinicalTrialTherapy;
DROP TABLE IF EXISTS clinicalTrialAgeGroup;
DROP TABLE IF EXISTS clinicalTrial;

DROP TABLE IF EXISTS drugGlobalApproavalStatusMolecularProfile;
DROP TABLE IF EXISTS drugGlobalApproavalStatusIndication;
DROP TABLE IF EXISTS drugGlobalApproavalStatusTherapy;
DROP TABLE IF EXISTS drugGlobalApproavalStatus;
DROP TABLE IF EXISTS drugTherapy;
DROP TABLE IF EXISTS drugEvidenceTreatmentApproch;
DROP TABLE IF EXISTS drugEvidenceReference;
DROP TABLE IF EXISTS drugEvidenceIndication;
DROP TABLE IF EXISTS drugEvidenceTherapy;
DROP TABLE IF EXISTS drugEvidenceMolecularProfile;
DROP TABLE IF EXISTS drugEvidence;
DROP TABLE IF EXISTS drugClinicalTrialTherapy;
DROP TABLE IF EXISTS drugClinicalTrial;
DROP TABLE IF EXISTS drugDrugClass;
DROP TABLE IF EXISTS drugDescriptionReference;
DROP TABLE IF EXISTS drugDescription;
DROP TABLE IF EXISTS drugSynonym;
DROP TABLE IF EXISTS drugTerm;
DROP TABLE IF EXISTS drug;

CREATE TABLE clinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    nctId varchar(255),
    title varchar(500),
    phase varchar(255),
    recruitment varchar(255),
    gender varchar(255),
    variantRequirement varchar(255),
    sponsors varchar(255),
    updateDate varchar(255),
    PRIMARY KEY (id)
);

CREATE TABLE clinicalTrialAgeGroup
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ageGroup varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

CREATE TABLE clinicalTrialTherapy
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    therapyId varchar(255),
    therapyName varchar(255),
    synonyms varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

CREATE TABLE clinicalTrialIndication
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    indicationId varchar(255),
    name varchar(255),
    source varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

CREATE TABLE clinicalTrialVariantRequirementDetail
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    requirementType varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

CREATE TABLE clinicalTrialVariantRequirementDetailMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialVariantRequirementDetailId int NOT NULL,
    idMolecularProfile varchar(255),
    profileName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialVariantRequirementDetailId) REFERENCES clinicalTrialVariantRequirementDetail(id)
);

CREATE TABLE clinicalTrialLocation
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    nctId varchar(255),
    facility varchar(255),
    city varchar(255),
    country varchar(255),
    status varchar(255),
    state varchar(255),
    zip varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

CREATE TABLE clinicalTrialContact
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialLocationId int NOT NULL,
    name varchar(255),
    email varchar(255),
    phone varchar(255),
    phoneExt varchar(255),
    role varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialLocationId) REFERENCES clinicalTrialLocation(id)
);

CREATE TABLE drug
(   id int NOT NULL AUTO_INCREMENT,
    idOfDrug varchar(255),
    drugName varchar(500),
    tradeName varchar(255),
    casRegistryNum varchar(255),
    nctId varchar(255),
    createDate varchar(255),
    PRIMARY KEY (id)
);

CREATE TABLE drugTerm
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    term varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugSynonym
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    synonym varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugDescription
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    description varchar(1500),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugDescriptionReference
(   id int NOT NULL AUTO_INCREMENT,
    drugDescriptionId int NOT NULL,
    drugDescriptionReferenceId varchar(255),
    pubMedId varchar(255),
    title varchar(500),
    url varchar(255),
    authors varchar(1000),
    journal varchar(255),
    volume varchar(255),
    issue varchar(255),
    date varchar(255),
    abstractText varchar(5000),
    year varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugDescriptionId) REFERENCES drugDescription(id)
);

CREATE TABLE drugDrugClass
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    drugDrugClassId varchar(255),
    drugClass varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugClinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    nctId varchar(255),
    title varchar(500),
    phase varchar(255),
    recruitment varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugClinicalTrialTherapy
(   id int NOT NULL AUTO_INCREMENT,
    drugClinicalTrialId int NOT NULL,
    drugClinicalTrialTherapyId varchar(255),
    therapyName varchar(255),
    synonym varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugClinicalTrialId) REFERENCES drugClinicalTrial(id)
);

CREATE TABLE drugEvidence
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    drugEvidenceId varchar(255),
    approvalStatus varchar(255),
    evidenceType varchar(255),
    efficacyEvidence varchar(1000),
    responseType varchar(255),
    ampCapAscoEvidenceLevel varchar(255),
    ampCapAscoInferredTier varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugEvidenceMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    drugEvidenceId int NOT NULL,
    drugEvidenceMolecularProfileId varchar(255),
    profileName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugEvidenceId) REFERENCES drugEvidence(id)
);

CREATE TABLE drugEvidenceTherapy
(   id int NOT NULL AUTO_INCREMENT,
    drugEvidenceId int NOT NULL,
    drugEvidenceTherapyId varchar(255),
    therapyName varchar(255),
    synonym varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugEvidenceId) REFERENCES drugEvidence(id)
);

CREATE TABLE drugEvidenceIndication
(   id int NOT NULL AUTO_INCREMENT,
    drugEvidenceId int NOT NULL,
    drugEvidenceIndicationId varchar(255),
    name varchar(255),
    source varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugEvidenceId) REFERENCES drugEvidence(id)
);

CREATE TABLE drugEvidenceReference
(   id int NOT NULL AUTO_INCREMENT,
    drugEvidenceId int NOT NULL,
    drugEvidenceReferenceId varchar(255),
    pubMedId varchar(255),
    title varchar(500),
    url varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugEvidenceId) REFERENCES drugEvidence(id)
);

CREATE TABLE drugEvidenceTreatmentApproch
(   id int NOT NULL AUTO_INCREMENT,
    drugEvidenceId int NOT NULL,
    drugEvidenceTreatmentApprochId varchar(255),
    name varchar(255),
    profileName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugEvidenceId) REFERENCES drugEvidence(id)
);

CREATE TABLE drugTherapy
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    drugTherapyId varchar(255),
    therapyname varchar(255),
    synonym varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugGlobalApproavalStatus
(   id int NOT NULL AUTO_INCREMENT,
    drugId int NOT NULL,
    drugGlobalApproavalStatusId varchar(255),
    approvalAuthority varchar(255),
    approvalStatus varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugId) REFERENCES drug(id)
);

CREATE TABLE drugGlobalApproavalStatusTherapy
(   id int NOT NULL AUTO_INCREMENT,
    drugGlobalApproavalStatusId int NOT NULL,
    drugGlobalApproavalStatusTherapyId varchar(255),
    therapyName varchar(255),
    synonym varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugGlobalApproavalStatusId) REFERENCES drugGlobalApproavalStatus(id)
);

CREATE TABLE drugGlobalApproavalStatusIndication
(   id int NOT NULL AUTO_INCREMENT,
    drugGlobalApproavalStatusId int NOT NULL,
    drugGlobalApproavalStatusIndicationId varchar(255),
    name varchar(255),
    source varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugGlobalApproavalStatusId) REFERENCES drugGlobalApproavalStatus(id)
);

CREATE TABLE drugGlobalApproavalStatusMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    drugGlobalApproavalStatusId int NOT NULL,
    drugGlobalApproavalStatusMolecularProfileId varchar(255),
    profileName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (drugGlobalApproavalStatusId) REFERENCES drugGlobalApproavalStatus(id)
);