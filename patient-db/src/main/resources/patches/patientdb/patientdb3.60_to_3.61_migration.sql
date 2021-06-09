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
    interpretation varchar(255) NOT NULL,
    reported BOOLEAN NOT NULL,
    PRIMARY KEY (id)
);