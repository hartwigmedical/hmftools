CREATE TABLE snomed
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    snomedId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);