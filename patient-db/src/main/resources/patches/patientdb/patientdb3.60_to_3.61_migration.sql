ALTER TABLE cuppa
    CHANGE COLUMN cuppaPrediction cuppaPrediction varchar(255);

CREATE TABLE virusAnnotation
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    taxid int NOT NULL,
    virusName varchar(255) NOT NULL,
    qcStatus varchar(255) NOT NULL,
    integrations int NOT NULL,
    interpretation varchar(255),
    reported BOOLEAN NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE copyNumberChromosomeArm
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    chromosomeArm varchar(255) NOT NULL,
    copyNumber DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);