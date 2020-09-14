CREATE TABLE amberAnonymous
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    sampleId varchar(255) NOT NULL,
    hmfSampleId varchar(255) NOT NULL,
    PRIMARY KEY (sampleId),
    INDEX(hmfSampleId)
);
