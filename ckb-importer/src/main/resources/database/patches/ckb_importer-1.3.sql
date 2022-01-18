## Adding relevant drug classes of therapies

DROP TABLE IF EXISTS treatmentApproachEvidence;
CREATE TABLE treatmentApproachEvidence
(   evidenceId int NOT NULL,
    treatmentApproachEvidenceId int NOT NULL,
    PRIMARY KEY (evidenceId, treatmentApproachEvidenceId),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id),
    FOREIGN KEY (treatmentApproachEvidenceId) REFERENCES treatmentApproach(id)
);

DROP TABLE IF EXISTS treatmentApproach;
CREATE TABLE treatmentApproach
(   id int NOT NULL AUTO_INCREMENT,
    treatmentApproachId int NOT NULL,
    createDate DATE NOT NULL,
    updateDate DATE,
    PRIMARY KEY (id)
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
