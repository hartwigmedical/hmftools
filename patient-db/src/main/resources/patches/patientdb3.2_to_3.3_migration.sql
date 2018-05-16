CREATE TABLE ranoMeasurement
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    responseDate DATE,
    therapyGiven varchar(50),
    targetLesionResponse varchar(50),
    noTargetLesionResponse varchar(50),
    overallResponse varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);