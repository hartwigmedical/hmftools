ALTER TABLE svFusion
    CHANGE skippedExons skippedExonsUp INT,
    ADD skippedExonsDown INT;

ALTER TABLE svFusion
    ADD fusedExonUp INT,
    ADD fusedExonDown INT;

ALTER TABLE svLink
    ADD ploidyUncertainty DOUBLE PRECISION,
    CHANGE chainIndex chainIndex VARCHAR(50);

CREATE TABLE svDriver
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NULL,
    gene VARCHAR(50) NOT NULL,
    eventType VARCHAR(50),
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(clusterId)
);