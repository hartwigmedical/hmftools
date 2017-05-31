CREATE TABLE patient
(   id int NOT NULL AUTO_INCREMENT,
    cpctId varchar(255) DEFAULT NULL,
    registrationDate DATE,
    gender varchar(10),
    ethnicity varchar(255),
    hospital varchar(255),
    birthYear int,
    primaryTumorLocation varchar(255),
    deathDate DATE,
    PRIMARY KEY (id)
);

CREATE TABLE sample
(   sampleId varchar(20) NOT NULL,
    patientId int NOT NULL,
    arrivalDate DATE NOT NULL,
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

CREATE TABLE biopsy
(   id int NOT NULL,
    sampleId varchar(20),
    patientId int NOT NULL,
    biopsyLocation varchar(255),
    biopsyDate DATE,
    PRIMARY KEY (id),
    FOREIGN KEY (sampleId) REFERENCES sample(sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

CREATE TABLE treatment
 (  id int NOT NULL,
    biopsyId int,
    patientId int NOT NULL,
    treatmentGiven varchar(3),
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (biopsyId) REFERENCES biopsy(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
 );

CREATE TABLE drug
 (  id int NOT NULL AUTO_INCREMENT,
    treatmentId int,
    patientId int NOT NULL,
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentId) REFERENCES treatment(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
 );

 CREATE TABLE treatmentResponse
  (  id int NOT NULL AUTO_INCREMENT,
     treatmentId int,
     patientId int NOT NULL,
     measurementDone varchar(5),
     responseDate DATE,
     response varchar(25),
     PRIMARY KEY (id),
     FOREIGN KEY (treatmentId) REFERENCES treatment(id),
     FOREIGN KEY (patientId) REFERENCES patient(id)
  );

CREATE TABLE somaticVariant
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(20) NOT NULL,
    patientId int NOT NULL,
    gene varchar(255) NOT NULL,
    position varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    cosmicId varchar(255),
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (sampleId) REFERENCES sample(sampleId),
    FOREIGN KEY (patientId) references patient(id)
);

CREATE TABLE comprehensiveSomaticVariant
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    position int not null,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

CREATE TABLE copyNumber
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    bafCount int not null,
    observedBaf DOUBLE PRECISION not null,
    actualBaf DOUBLE PRECISION not null,
    copyNumber DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

CREATE TABLE purity
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    purity DOUBLE PRECISION not null,
    normFactor DOUBLE PRECISION not null,
    score DOUBLE PRECISION not null,
    ploidy DOUBLE PRECISION not null,
    diploidProportion DOUBLE PRECISION not null,
    polyclonalProportion DOUBLE PRECISION not null,
    minPurity DOUBLE PRECISION not null,
    maxPurity DOUBLE PRECISION not null,
    minPloidy DOUBLE PRECISION not null,
    maxPloidy DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);
