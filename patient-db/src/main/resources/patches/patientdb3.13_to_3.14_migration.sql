CREATE TABLE clinicalEvidence
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(255) NOT NULL,
    eventType varchar(255) NOT NULL,
    eventMatch varchar(255) NOT NULL,
    nameEvidence varchar(500) NOT NULL,
    typeEvidence varchar(255) NOT NULL,
    response varchar(255) NOT NULL,
    levelEvidence varchar(50) NOT NULL,
    sourceEvidence varchar(255) NOT NULL,
    cancerType varchar(500) NOT NULL,
    isOnLabel BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);