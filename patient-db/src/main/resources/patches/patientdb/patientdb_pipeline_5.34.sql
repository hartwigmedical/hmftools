####
# SQL updates for Pipeline release 5.34
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# ACTIN-341

DROP TABLE IF EXISTS svBreakend;
CREATE TABLE svBreakend
(   id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId VARCHAR(50) NOT NULL,
    svId INT NOT NULL,
    startBreakend TINYINT(1) NOT NULL,
    gene VARCHAR(50) NOT NULL, -- # length here comes from ensembl db schema
    transcriptId VARCHAR(128) NOT NULL, -- # length here comes from ensembl db schema
    canonicalTranscript TINYINT(1) NOT NULL,
    geneOrientation VARCHAR(20) NOT NULL,
    disruptive TINYINT(1) NOT NULL,
    reportedDisruption TINYINT(1) NOT NULL,
    undisruptedCopyNumber DOUBLE PRECISION,
    regionType VARCHAR(20) NOT NULL,
    codingContext VARCHAR(20),
    biotype VARCHAR(255),
    exonUp SMALLINT NOT NULL,
    exonDown SMALLINT NOT NULL,
    exonicBasePhase TINYINT,
    nextSpliceExonRank SMALLINT,
    nextSpliceExonPhase TINYINT,
    nextSpliceDistance INT,
    totalExonCount SMALLINT NOT NULL,
    PRIMARY KEY (id)
);

CREATE INDEX svBreakend_sampleId_svId ON svBreakend (sampleId, svId);
CREATE INDEX svBreakend_gene ON svBreakend (gene);
CREATE INDEX svBreakend_transcriptId ON svBreakend (transcriptId);
