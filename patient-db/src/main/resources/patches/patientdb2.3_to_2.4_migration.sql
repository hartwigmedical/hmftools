DROP TABLE IF EXISTS patient;
CREATE TABLE patient
(   id int NOT NULL AUTO_INCREMENT,
    patientIdentifier varchar(50),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS baseline;
CREATE TABLE baseline
(   patientId int NOT NULL,
    registrationDate DATE,
    informedConsentDate DATE,
    gender varchar(10),
    hospital varchar(255),
    birthYear int,
    cancerType varchar(255),
    cancerSubtype varchar(255),
    deathDate DATE,
    hasSystemicPreTreatment varchar(3),
    hasRadiotherapyPreTreatment varchar(3),
    preTreatments varchar(255),
    PRIMARY KEY (patientId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);