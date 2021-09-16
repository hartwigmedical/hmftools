DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;