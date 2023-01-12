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
    reportGermlineDisruption BOOLEAN NOT NULL,
    additionalReportedTranscripts varchar(255) NOT NULL,
    reportPGX BOOLEAN NOT NULL,
    PRIMARY KEY (gene)
);

# Germline SVs

ALTER TABLE structuralVariantGermline
    DROP COLUMN gene,
    ADD COLUMN homologySequenceStart varchar(255) not null after qualScore,
    ADD COLUMN homologySequenceEnd varchar(255) after startHomologySequence,
    ADD COLUMN junctionCopyNumber DOUBLE PRECISION after endHomologySequence,
    ADD COLUMN adjustedAFStart DOUBLE PRECISION after junctionCopyNumber,
    ADD COLUMN adjustedAFEnd DOUBLE PRECISION after adjustedAFStart,
    ADD COLUMN adjustedCopyNumberStart DOUBLE PRECISION after adjustedAFEnd,
    ADD COLUMN adjustedCopyNumberEnd DOUBLE PRECISION after adjustedCopyNumberStart,
    ADD COLUMN adjustedCopyNumberChangeStart DOUBLE PRECISION after adjustedCopyNumberEnd,
    ADD COLUMN adjustedCopyNumberChangeEnd DOUBLE PRECISION after adjustedCopyNumberChangeStart;

# Isofox
ALTER TABLE rnaStatistics
	CHANGE COLUMN qcStatus qcStatus varchar(100);


