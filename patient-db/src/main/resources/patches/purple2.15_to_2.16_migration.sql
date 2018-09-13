DROP TABLE IF EXISTS driverCatalog;
CREATE TABLE driverCatalog
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    category varchar(255) NOT NULL,
    driver varchar(255) NOT NULL,
    driverLikelihood DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(gene)
);
