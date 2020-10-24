ALTER TABLE baseline
    DROP COLUMN cancerSubtype,
    ADD COLUMN primaryTumorSubLocation varchar(255) AFTER primaryTumorLocation,
    ADD COLUMN primaryTumorType varchar(255) AFTER primaryTumorSubLocation,
    ADD COLUMN primaryTumorSubType varchar(255) AFTER primaryTumorType,
    ADD COLUMN primaryTumorExtraDetails varchar(255) AFTER primaryTumorSubType,
    ADD COLUMN primaryTumorOverridden BOOLEAN AFTER primaryTumorExtraDetails

CREATE TABLE doidEntry
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    doid varchar(255) NOT NULL,
    doidTerm varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);