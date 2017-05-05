CREATE TABLE patientInfo
(   id int NOT NULL AUTO_INCREMENT,
    cpctId varchar(255) DEFAULT NULL,
    registrationDate DATE,
    gender varchar(10),
    ethnicity varchar(255),
    hospital varchar(255),
    birthYear int,
    tumorLocation varchar(255),
    deathDate DATE,
    PRIMARY KEY (id))
);


CREATE TABLE biopsyLimsData
(   sampleId varchar(20) NOT NULL,
    arrivalDate DATE NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patientInfo(id)
);

CREATE TABLE biopsyClinicalData
(   id int NOT NULL,
    location varchar(255),
    date DATE,
    sampleId varchar(20),
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patientInfo(id),
    FOREIGN KEY (sampleId) REFERENCES biopsyLimsData(sampleId)
);

CREATE TABLE treatmentData
 (  id int NOT NULL,
    treatmentGiven varchar(3),
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    clinicalBiopsyId int,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patientInfo(id),
    FOREIGN KEY (clinicalBiopsyId) REFERENCES biopsyClinicalData(id)
 );

CREATE TABLE drugData
 (  id int NOT NULL AUTO_INCREMENT,
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    treatmentId int,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patientInfo(id),
    FOREIGN KEY (treatmentId) REFERENCES treatmentData(id)
 );

 CREATE TABLE treatmentResponseData
  (  id int NOT NULL AUTO_INCREMENT,
     date DATE,
     response varchar(25),
     measurementDone varchar(5),
     treatmentId int,
     patientId int NOT NULL,
     PRIMARY KEY (id),
     FOREIGN KEY (patientId) REFERENCES patientInfo(id),
     FOREIGN KEY (treatmentId) REFERENCES treatmentData(id)
  );

CREATE TABLE somaticVariantData
(   id int NOT NULL AUTO_INCREMENT,
    gene varchar(255) NOT NULL,
    position varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    cosmicId varchar(255),
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) references patientInfo(id)
);
