ALTER TABLE driverCatalog
    ADD biallelic BOOLEAN NOT NULL AFTER inframe;

ALTER TABLE driverCatalog
    CHANGE driver likelihoodMethod varchar(255) NOT NULL,
    ADD driver varchar(255) NOT NULL AFTER category;

UPDATE driverCatalog set driver = "MUTATION";
UPDATE driverCatalog set driver = "AMP" WHERE likelihoodMethod = "AMP";
UPDATE driverCatalog set driver = "DEL" WHERE likelihoodMethod = "DEL";

ALTER TABLE driverCatalog DROP INDEX sampleId;
ALTER TABLE driverCatalog ADD INDEX (sampleId, gene);

ALTER TABLE geneCopyNumber DROP INDEX sampleId;
ALTER TABLE geneCopyNumber ADD INDEX (sampleId, gene);