DROP TABLE IF EXISTS driverCatalog;
CREATE TABLE driverCatalog
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    category varchar(255) NOT NULL,
    driver varchar(255) NOT NULL,
    dndsLikelihood DOUBLE PRECISION NOT NULL,
    driverLikelihood DOUBLE PRECISION NOT NULL,
    missense int NOT NULL,
    nonsense int NOT NULL,
    splice int NOT NULL,
    frameshift int NOT NULL,
    inframe int NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(gene)
);