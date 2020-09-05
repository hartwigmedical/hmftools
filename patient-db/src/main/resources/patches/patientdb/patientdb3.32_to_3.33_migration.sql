DROP TABLE IF EXISTS rna;
CREATE TABLE rna
(   sampleId varchar(255) NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS ecrf;
DROP TABLE IF EXISTS ecrfDatamodel;
SET FOREIGN_KEY_CHECKS = 0;
DROP TABLE IF EXISTS sample;
SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE cpctEcrf
(   id int NOT NULL AUTO_INCREMENT,
    patientId varchar(20),
    studyEvent varchar(100),
    studyEventKey int not null,
    form varchar(100),
    formKey int not null,
    itemGroup varchar(100),
    itemGroupKey int not null,
    item varchar(100),
    itemValue varchar(1500),
    status varchar(30),
    locked varchar(5),
    sequenced varchar(5),
    fieldName varchar(100),
    relevant varchar(5),
    PRIMARY KEY (id),
    INDEX(patientId),
    INDEX(studyEvent),
    INDEX(form),
    INDEX(itemGroup),
    INDEX(item),
    INDEX(itemValue (255)),
    INDEX(status),
    INDEX(locked),
    INDEX(sequenced),
    INDEX(fieldName),
    INDEX(relevant)
);

CREATE TABLE cpctEcrfDatamodel
(   fieldName varchar(100),
    description varchar(500),
    codeList varchar(3000),
    relevant varchar(5)
);

CREATE TABLE sample
(   sampleId varchar(255) NOT NULL,
    patientId int NOT NULL,
    setName varchar(255) NOT NULL,
    arrivalDate DATE NOT NULL,
    samplingDate DATE,
    dnaNanograms int,
    limsPrimaryTumor varchar(255),
    pathologyTumorPercentage varchar(100),
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);
