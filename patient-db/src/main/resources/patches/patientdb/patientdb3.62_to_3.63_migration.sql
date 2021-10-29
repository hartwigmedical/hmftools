DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;

ALTER TABLE baseline
    ADD COLUMN pifVersion varchar(255) AFTER informedConsentDate;

ALTER TABLE baseline
    ADD COLUMN inHMFDatabase BOOLEAN AFTER pifVersion;

ALTER TABLE baseline
    ADD COLUMN outsideEU BOOLEAN AFTER inHMFDatabase;

ALTER TABLE virusAnnotation
    ADD COLUMN percentageCovered BOOLEAN AFTER interpretation;

ALTER TABLE virusAnnotation
    ADD COLUMN meanCoverage BOOLEAN AFTER percentageCovered;

ALTER TABLE virusAnnotation
    ADD COLUMN expectedClonalCoverage BOOLEAN AFTER meanCoverage;
