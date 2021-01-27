SET FOREIGN_KEY_CHECKS = 0;

-- TODO Remove per 1st of March 2021
DROP TABLE IF EXISTS ckbEntry;

DROP TABLE IF EXISTS clinicalTrial;
CREATE TABLE clinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    nctId varchar(255),
    title varchar(255),
    phase varchar(255),
    recruitment varchar(255),
    gender varchar(255),
    variantRequirement varchar(255),
    sponsors varchar(255),
    updateDate varchar(255),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS clinicalTrialAgeGroup;
CREATE TABLE clinicalTrialAgeGroup
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    ageGroup varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialTherapy;
CREATE TABLE clinicalTrialTherapy
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    therapyId varchar(255),
    therapyName varchar(255),
    synonyms varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialIndication;
CREATE TABLE clinicalTrialIndication
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    indicationId varchar(255),
    name varchar(255),
    source varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetail;
CREATE TABLE clinicalTrialVariantRequirementDetail
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialId int NOT NULL,
    requirementType varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialId) REFERENCES clinicalTrial(id)
);

DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetailMolecularProfile;
CREATE TABLE clinicalTrialVariantRequirementDetailMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    clinicalTrialVariantRequirementDetailId int NOT NULL,
    idMolecularProfile varchar(255),
    profileName varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (clinicalTrialVariantRequirementDetailId) REFERENCES clinicalTrialVariantRequirementDetail(id)
);

DROP TABLE IF EXISTS clinicalTrialLocation;
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

DROP TABLE IF EXISTS clinicalTrialContact;
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