DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;

ALTER TABLE virusAnnotation
    ADD COLUMN percentageCovered double precision not null AFTER interpretation;

ALTER TABLE virusAnnotation
    ADD COLUMN meanCoverage double precision not null AFTER percentageCovered;

ALTER TABLE virusAnnotation
    ADD COLUMN expectedClonalCoverage double precision AFTER meanCoverage;