DROP TABLE if EXISTS bqr;
CREATE TABLE bqr
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(255) NOT NULL,
    sampleType varchar(50) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    trinucleotideContext varchar(3) NOT NULL,
    bqrCount int NOT NULL,
    origQuality DOUBLE PRECISION NOT NULL,
    recalibratedQuality DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id)
);