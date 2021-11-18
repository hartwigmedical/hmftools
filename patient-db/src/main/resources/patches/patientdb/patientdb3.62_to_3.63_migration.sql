DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;

ALTER TABLE baseline
    ADD COLUMN pifVersion varchar(255) AFTER informedConsentDate;

ALTER TABLE baseline
    ADD COLUMN inHMFDatabase BOOLEAN AFTER pifVersion;

ALTER TABLE baseline
    ADD COLUMN outsideEU BOOLEAN AFTER inHMFDatabase;

ALTER TABLE copyNumber
	ADD UNIQUE KEY (sampleId, chromosome, start, end);

ALTER TABLE geneCopyNumber
	MODIFY COLUMN sampleId VARCHAR(31);
ALTER TABLE geneCopyNumber
	MODIFY COLUMN geneCopyNumber VARCHAR(31);
ALTER TABLE geneCopyNumber
	ADD UNIQUE KEY (sampleId, chromosome, gene, transcriptId);

ALTER TABLE somaticVariant
	MODIFY COLUMN sampleId VARCHAR(31);
ALTER TABLE somaticVariant
	MODIFY COLUMN chromosome VARCHAR(31);
ALTER TABLE somaticVariant
	ADD UNIQUE KEY (sampleId, chromosome, position, ref, alt);

ALTER TABLE virusAnnotation
    ADD COLUMN percentageCovered double NOT NULL AFTER interpretation;

ALTER TABLE virusAnnotation
    ADD COLUMN meanCoverage double NOT NULL AFTER percentageCovered;

ALTER TABLE virusAnnotation
    ADD COLUMN expectedClonalCoverage double AFTER meanCoverage;