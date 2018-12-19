CREATE TABLE clinicalTrials
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(255) NOT NULL,
    eventType varchar(255) NOT NULL,
    eventMatch varchar(255) NOT NULL,
    trial varchar(500) NOT NULL,
    cancerType varchar(500) NOT NULL,
    label varchar(50) NOT NULL,
    evidenceSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

CREATE TABLE clinicalEvidence
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(255) NOT NULL,
    eventType varchar(255) NOT NULL,
    eventMatch varchar(255) NOT NULL,
    drug varchar(500) NOT NULL,
    drugsType varchar(255) NOT NULL,
    response varchar(255) NOT NULL,
    cancerType varchar(500) NOT NULL,
    label varchar(50) NOT NULL,
    evidenceLevel varchar(50) NOT NULL,
    evidenceSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);