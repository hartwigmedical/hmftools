CREATE TABLE snomed
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    snomedConceptId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS doidNode;
CREATE TABLE doidNode
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    doid varchar(255) NOT NULL,
    doidTerm varchar(255) NOT NULL,
    snomedConceptId varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);
