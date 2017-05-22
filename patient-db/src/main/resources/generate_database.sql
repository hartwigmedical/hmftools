CREATE TABLE patients
(   id int NOT NULL AUTO_INCREMENT,
    cpctId varchar(255) DEFAULT NULL,
    registrationDate DATE,
    gender varchar(10),
    ethnicity varchar(255),
    hospital varchar(255),
    birthYear int,
    tumorLocation varchar(255),
    deathDate DATE,
    PRIMARY KEY (id)
);


CREATE TABLE limsBiopsies
(   sampleId varchar(20) NOT NULL,
    arrivalDate DATE NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patients(id)
);

CREATE TABLE clinicalBiopsies
(   id int NOT NULL,
    location varchar(255),
    date DATE,
    sampleId varchar(20),
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patients(id),
    FOREIGN KEY (sampleId) REFERENCES limsBiopsies(sampleId)
);

CREATE TABLE treatments
 (  id int NOT NULL,
    treatmentGiven varchar(3),
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    clinicalBiopsyId int,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patients(id),
    FOREIGN KEY (clinicalBiopsyId) REFERENCES clinicalBiopsies(id)
 );

CREATE TABLE drugs
 (  id int NOT NULL AUTO_INCREMENT,
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    treatmentId int,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patients(id),
    FOREIGN KEY (treatmentId) REFERENCES treatments(id)
 );

 CREATE TABLE treatmentResponses
  (  id int NOT NULL AUTO_INCREMENT,
     date DATE,
     response varchar(25),
     measurementDone varchar(5),
     treatmentId int,
     patientId int NOT NULL,
     PRIMARY KEY (id),
     FOREIGN KEY (patientId) REFERENCES patients(id),
     FOREIGN KEY (treatmentId) REFERENCES treatments(id)
  );

CREATE TABLE somaticVariants
(   id int NOT NULL AUTO_INCREMENT,
    gene varchar(255) NOT NULL,
    position varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    cosmicId varchar(255),
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    sampleId varchar(20) NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) references patients(id),
    FOREIGN KEY (sampleId) REFERENCES limsBiopsies(sampleId)
);
