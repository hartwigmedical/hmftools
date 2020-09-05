ALTER TABLE structuralVariant
    ADD COLUMN svId INT NOT NULL AFTER modified,
    CHANGE adjustedStartAF adjustedAFStart DOUBLE PRECISION,
    CHANGE adjustedEndAF adjustedAFEnd DOUBLE PRECISION,
    CHANGE adjustedStartCopyNumber adjustedCopyNumberStart DOUBLE PRECISION,
    CHANGE adjustedEndCopyNumber adjustedCopyNumberEnd DOUBLE PRECISION,
    CHANGE adjustedStartCopyNumberChange adjustedCopyNumberChangeStart DOUBLE PRECISION,
    CHANGE adjustedEndCopyNumberChange adjustedCopyNumberChangeEnd DOUBLE PRECISION,
    ADD INDEX (svId);

UPDATE structuralVariant
SET svId = id;

CREATE TABLE svAnnotation
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    svId INT NOT NULL,
    clusterId INT NOT NULL,
    clusterReason VARCHAR(255) NULL,
    fragileSiteStart BOOLEAN NOT NULL,
    fragileSiteEnd BOOLEAN NOT NULL,
    isFoldback BOOLEAN NOT NULL,
    lineTypeStart VARCHAR(20),
    lineTypeEnd VARCHAR(20),
    ploidyMin DOUBLE PRECISION,
    ploidyMax DOUBLE PRECISION,
    geneStart VARCHAR(100),
    geneEnd varchar(100),
    replicationTimingStart DOUBLE PRECISION,
    replicationTimingEnd DOUBLE PRECISION,
    localTopologyIdStart INT,
    localTopologyIdEnd INT,
    localTopologyStart varchar(20),
    localTopologyEnd VARCHAR(20),
    localTICountStart INT,
    localTICountEnd INT,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(clusterId),
    INDEX(svId)
);

CREATE TABLE svCluster
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NOT NULL,
    resolvedType VARCHAR(20),
    synthetic BOOLEAN NOT NULL,
    subClonal BOOLEAN NOT NULL,
    subType VARCHAR(20),
    clusterCount INT,
    clusterDesc VARCHAR(50),
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(clusterId)
);

CREATE TABLE svLink
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NOT NULL,
    chainId INT NOT NULL,
    chainIndex INT NOT NULL,
    chainLinkCount INT NOT NULL,
    lowerSvId INT NOT NULL,
    upperSvId INT NOT NULL,
    lowerBreakendIsStart BOOLEAN NOT NULL,
    upperBreakendIsStart BOOLEAN NOT NULL,
    chromosome VARCHAR(3),
    arm VARCHAR(3),
    assembled BOOLEAN NOT NULL,
    traversedSVCount INT,
    linkLength INT,
    ploidy DOUBLE PRECISION,
    pseudogeneInfo varchar(255),
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(clusterId)
);

RENAME TABLE structuralVariantBreakend TO svBreakend;

ALTER TABLE svBreakend
    CHANGE isStartEnd startBreakend BOOLEAN,
    CHANGE isCanonicalTranscript canonicalTranscript BOOLEAN,
    CHANGE structuralVariantId svId INT,
    DROP COLUMN geneId,
    DROP COLUMN strand,
    DROP COLUMN exonRankUpstream,
    DROP COLUMN exonRankDownstream,
    DROP COLUMN exonPhaseUpstream,
    DROP COLUMN exonPhaseDownstream,
    ADD geneOrientation VARCHAR(20) NOT NULL,
    ADD disruptive BOOLEAN NOT NULL,
    ADD reportedDisruption BOOLEAN NOT NULL,
    ADD regionType VARCHAR(20) NOT NULL,
    ADD codingContext VARCHAR(20),
    ADD biotype VARCHAR(255),
    ADD exonicBasePhase TINYINT,
    ADD nextSpliceExonRank TINYINT UNSIGNED,
    ADD nextSpliceExonPhase TINYINT,
    ADD nextSpliceDistance INT,
    CHANGE exonMax totalExonCount SMALLINT;

RENAME TABLE structuralVariantFusion TO svFusion;

ALTER TABLE svFusion
	ADD name VARCHAR(50) NOT NULL,
    CHANGE isReported reported BOOLEAN,
    ADD reportedType varchar(255) NULL,
    ADD chainLength INT,
    ADD chainLinks INT,
    ADD chainTerminated BOOLEAN,
    ADD domainsKept VARCHAR(255),
    ADD domainsLost VARCHAR(255),
    ADD skippedExons INT;

DROP TABLE structuralVariantDisruption;

CREATE TABLE viralInsertion
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    svId INT NOT NULL,
    virusId VARCHAR(50) NOT NULL,
    virusName VARCHAR(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(svId)
);
