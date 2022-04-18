####
# SQL updates for Pipeline release 5.28 which only support PROTECT v2.1
# PROTECT

DROP TABLE IF EXISTS protect;
CREATE TABLE protect
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255),
    transcript varchar(255),
    isCanonical BOOLEAN,
    event varchar(255) NOT NULL,
    eventIsHighDriver BOOLEAN,
    germline BOOLEAN NOT NULL,
    reported BOOLEAN NOT NULL,
    treatment varchar(255) NOT NULL,
    onLabel BOOLEAN NOT NULL,
    level varchar(255) NOT NULL,
    direction varchar(255) NOT NULL,
    source varchar(255) NOT NULL,
    sourceEvent varchar(255) NOT NULL,
    sourceUrls varchar(2500) NOT NULL,
    evidenceType varchar(50) NOT NULL,
    rangeRank int,
    evidenceUrls varchar(2500) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);