DROP TABLE IF EXISTS amberPatient;
CREATE TABLE amberPatient
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    sampleId varchar(255) NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS amberMapping;
CREATE TABLE amberMapping
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    firstSampleId varchar(255) NOT NULL,
    secondSampleId varchar(255) NOT NULL,
    matches int NOT NULL,
    sites int NOT NULL,
    likelihood DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (firstSampleId, secondSampleId)
);

DROP TABLE IF EXISTS amber;
