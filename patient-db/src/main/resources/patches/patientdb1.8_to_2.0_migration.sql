ALTER TABLE patient
    ADD COLUMN informedConsentDate DATE AFTER registrationDate;

ALTER TABLE patient
    ADD COLUMN hasSystemicPreTreatment varchar(3) AFTER deathDate;

ALTER TABLE patient
    ADD COLUMN hasRadiotherapyPreTreatment varchar(3) AFTER hasSystemicPreTreatment;

ALTER TABLE patient
    ADD COLUMN preTreatments varchar(255) AFTER hasRadiotherapyPreTreatment;

CREATE TABLE preTreatmentDrug
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    bestResponse varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);