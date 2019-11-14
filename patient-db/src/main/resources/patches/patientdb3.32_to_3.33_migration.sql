DROP TABLE IF EXISTS rna;
CREATE TABLE rna
(   sampleId varchar(255) NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS cpctEcrf;
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

DROP TABLE IF EXISTS cpctEcrfDatamodel;
CREATE TABLE cpctEcrfDatamodel
(   fieldName varchar(100),
    description varchar(500),
    codeList varchar(3000),
    relevant varchar(5)
);