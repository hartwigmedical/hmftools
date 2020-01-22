DROP TABLE IF EXISTS clinicalEvidenceProtect;
CREATE TABLE clinicalEvidenceProtect
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    event varchar(255) NOT NULL,
    name varchar(500) NOT NULL,
    type varchar(255) NOT NULL,
    response varchar(255) NOT NULL,
    level varchar(50) NOT NULL,
    source varchar(255) NOT NULL,
    cancerType varchar(500) NOT NULL,
    isOnLabel BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);