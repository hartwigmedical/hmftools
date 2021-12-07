####
# RNA and Isofox
- rna table no longer used, replaced by rnaStatistics

DROP TABLE rna;

ALTER TABLE rnaStatistics
    ADD COLUMN qcStatus varchar(50) NOT NULL AFTER sampleId;

####
# Baseline changes
ALTER TABLE baseline
    ADD COLUMN pifVersion varchar(255) AFTER informedConsentDate;

ALTER TABLE baseline
    ADD COLUMN inHMFDatabase BOOLEAN AFTER pifVersion;

ALTER TABLE baseline
    ADD COLUMN outsideEU BOOLEAN AFTER inHMFDatabase;

####
# Added unique composite keys to some tables to avoid duplications
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

####
# Linx v1.17
# - removal of old viral insertion table
# - remove unused replication timing fields
# - add new germlineSV table

DROP TABLE viralInsertion;

ALTER TABLE svAnnotation
    DROP COLUMN replicationTimingStart,
    DROP COLUMN replicationTimingEnd;

CREATE TABLE structuralVariantGermline
(   id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosomeStart varchar(10) NOT NULL,
    chromosomeEnd varchar(10),
    positionStart int not null,
    positionEnd int,
    orientationStart tinyint not null,
    orientationEnd tinyint,
    gene varchar(50) NOT NULL,
    type varchar(10) NOT NULL,
    filter varchar(50) NOT NULL,
    event varchar(50),
    qualScore DOUBLE PRECISION,
    germlineFragments int,
    germlineReferenceFragmentsStart int,
    germlineReferenceFragmentsEnd int,
    tumorFragments int,
    tumorReferenceFragmentsStart int,
    tumorReferenceFragmentsEnd int,
    insertSequence varchar(2048) not null,
    insertSequenceAlignments varchar(512),
    insertSequenceRepeatClass varchar(64),
    insertSequenceRepeatType varchar(64),
    clusterId int null,
    clusterCount int null,
    resolvedType VARCHAR(20),
    linkedByStart varchar(1024),
    linkedByEnd varchar(1024),
    cohortFrequency int not null,
    reported BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);


####
# Purple v3.2
# - introduce germline deletion table and drivers
# - support alternative (ie non-canonical) transcripts for drivers and somatic variants
# - clean-up of impact fields with switch from SnpEff to Pave

ALTER TABLE geneCopyNumber
    ADD COLUMN canonicalTranscript BOOLEAN NOT NULL DEFAULT 1,
    DROP COLUMN germlineHomDeletionRegions,
    DROP COLUMN germlineHetToHomDeletionRegions;

DROP TABLE IF EXISTS copyNumberGermline;

ALTER TABLE somaticVariant
    ADD COLUMN spliceRegion BOOLEAN NOT NULL DEFAULT 0 AFTER canonicalHgvsProteinImpact,
    ADD COLUMN otherTranscriptEffects varchar(255) AFTER spliceRegion,
    CHANGE COLUMN genesEffected genesAffected int not null,
    DROP COLUMN worstEffect,
    DROP COLUMN worstEffectTranscript;

ALTER TABLE driverCatalog
    ADD COLUMN transcriptId varchar(50) NOT NULL AFTER gene,
    ADD COLUMN canonicalTranscript BOOLEAN NOT NULL DEFAULT 1 after transcriptId;

ALTER TABLE germlineVariant
    ADD COLUMN spliceRegion BOOLEAN NOT NULL DEFAULT 0 AFTER canonicalHgvsProteinImpact,
    ADD COLUMN otherTranscriptEffects varchar(255) AFTER spliceRegion,
    CHANGE COLUMN genesEffected genesAffected int not null,
    DROP COLUMN worstEffect,
    DROP COLUMN worstEffectTranscript;

CREATE TABLE germlineDeletion
(   id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(50) NOT NULL,
    chromosome varchar(10) NOT NULL,
    regionStart int not null,
    regionEnd int not null,
    depthWindowCount int not null,
    exonStart int not null,
    exonEnd int not null,
    detectionMethod varchar(20) NOT NULL,
    germlineStatus varchar(20) NOT NULL,
    tumorStatus varchar(20) NOT NULL,
    germlineCopyNumber DOUBLE PRECISION NOT NULL,
    tumorCopyNumber DOUBLE PRECISION NOT NULL,
    filter varchar(50) NOT NULL,
    cohortFrequency int not null,
    reported BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);

