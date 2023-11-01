####
# SQL updates for Pipeline release 5.32
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Driver gene panel

DROP TABLE IF EXISTS driverGenePanel;
CREATE TABLE driverGenePanel
(   modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    reportMissense BOOLEAN NOT NULL,
    reportNonsense BOOLEAN NOT NULL,
    reportSplice BOOLEAN NOT NULL,
    reportDeletion BOOLEAN NOT NULL,
    reportDisruption BOOLEAN NOT NULL,
    reportAmplification BOOLEAN NOT NULL,
    reportSomaticHotspot BOOLEAN NOT NULL,
    reportGermlineVariant varchar(50) NOT NULL,
    reportGermlineHotspot varchar(50) NOT NULL,
    likelihoodType varchar(255) NOT NULL,
    reportGermlineDisruption varchar(50) NOT NULL,
    reportGermlineDeletion varchar(50) NOT NULL,
    additionalReportedTranscripts varchar(255) NOT NULL,
    reportPGX BOOLEAN NOT NULL,
    PRIMARY KEY (gene)
);

# Germline SVs

ALTER TABLE structuralVariantGermline
    DROP COLUMN gene,
    ADD COLUMN svId INT NOT NULL after sampleId,
    ADD COLUMN vcfId varchar(50) NOT NULL after svId,
    ADD COLUMN homologySequenceStart varchar(255) not null after qualScore,
    ADD COLUMN homologySequenceEnd varchar(255) after homologySequenceStart,
    ADD COLUMN junctionCopyNumber DOUBLE PRECISION after homologySequenceEnd,
    ADD COLUMN adjustedAFStart DOUBLE PRECISION after junctionCopyNumber,
    ADD COLUMN adjustedAFEnd DOUBLE PRECISION after adjustedAFStart,
    ADD COLUMN adjustedCopyNumberStart DOUBLE PRECISION after adjustedAFEnd,
    ADD COLUMN adjustedCopyNumberEnd DOUBLE PRECISION after adjustedCopyNumberStart,
    ADD COLUMN adjustedCopyNumberChangeStart DOUBLE PRECISION after adjustedCopyNumberEnd,
    ADD COLUMN adjustedCopyNumberChangeEnd DOUBLE PRECISION after adjustedCopyNumberChangeStart,
    DROP COLUMN reported;

# Germline Breakends

CREATE TABLE svBreakendGermline
(   id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(50) NOT NULL,
    svId INT NOT NULL,
    startBreakend BOOLEAN NOT NULL,
    gene VARCHAR(50) NOT NULL, # length here comes from ensembl db schema
    transcriptId VARCHAR(128) NOT NULL, # length here comes from ensembl db schema
    canonicalTranscript BOOLEAN NOT NULL,
    geneOrientation VARCHAR(20) NOT NULL,
    disruptive BOOLEAN NOT NULL,
    reportedDisruption BOOLEAN NOT NULL,
    undisruptedCopyNumber DOUBLE PRECISION,
    regionType VARCHAR(20) NOT NULL,
    codingType VARCHAR(20),
    biotype VARCHAR(255),
    exonicBasePhase TINYINT,
    nextSpliceExonRank TINYINT UNSIGNED,
    nextSpliceExonPhase TINYINT,
    nextSpliceDistance INT,
    totalExonCount SMALLINT NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, svId),
    INDEX(gene),
    INDEX(transcriptId)
);

# Isofox
ALTER TABLE rnaStatistics
	CHANGE COLUMN qcStatus qcStatus varchar(100);


