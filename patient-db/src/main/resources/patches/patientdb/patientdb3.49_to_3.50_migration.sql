DROP TABLE IF EXISTS doidEntry;

CREATE TABLE doidNode
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    doid varchar(255) NOT NULL,
    doidTerm varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);