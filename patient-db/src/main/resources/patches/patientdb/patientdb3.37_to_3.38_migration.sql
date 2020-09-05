DROP TABLE IF EXISTS pgxCalls;
CREATE TABLE pgxCalls
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    positionGRCh37 varchar(255) NOT NULL,
    refGRCh37 varchar(255) NOT NULL,
    altGRCh37 varchar(255) NOT NULL,
    positionGRCh38 varchar(255) NOT NULL,
    refGRCh38 varchar(255) NOT NULL,
    altGRCh38 varchar(255) NOT NULL,
    rsid varchar(255) NOT NULL,
    variantAnnotation varchar(255) NOT NULL,
    filter varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS pgxGenotype;
CREATE TABLE pgxGenotype
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    haplotype varchar(255) NOT NULL,
    function varchar(255) NOT NULL,
    linkedDrugs varchar(255) NOT NULL,
    urlPrescriptionInfo varchar(255) NOT NULL,
    panelVersion varchar(255) NOT NULL,
    repoVersion varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

ALTER TABLE metric
    ADD metricQC BOOLEAN NOT NULL after tumorCoverage60xPercentage;