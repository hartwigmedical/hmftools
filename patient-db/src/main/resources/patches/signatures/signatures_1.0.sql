CREATE TABLE signature
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    signature VARCHAR(50) NOT NULL,
    allocation DOUBLE PRECISION,
    percent DOUBLE PRECISION,
    PRIMARY KEY (id),
    INDEX(sampleId)
);
