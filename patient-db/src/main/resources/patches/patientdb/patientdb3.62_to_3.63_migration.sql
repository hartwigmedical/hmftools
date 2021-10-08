DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;

ALTER TABLE baseline
    ADD COLUMN informedConsentDate varchar(255) AFTER pifVersion;

ALTER TABLE baseline
    ADD COLUMN pifVersion BOOLEAN AFTER inHMFDatabase;

ALTER TABLE baseline
    ADD COLUMN inHMFDatabase BOOLEAN AFTER outsideEU;
