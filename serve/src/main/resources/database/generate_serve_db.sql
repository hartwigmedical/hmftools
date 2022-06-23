SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS actionableHotspots;
CREATE TABLE actionableHotspots
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    chromosome varchar(50) NOT NULL,
    position int NOT NULL,
    ref varchar(100) NOT NULL,
    alt varchar(100) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableRanges;
CREATE TABLE actionableRanges
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    transcript varchar(50) NOT NULL,
    chromosome varchar(50) NOT NULL,
    start int NOT NULL,
    end int NOT NULL,
    mutationType varchar(50) NOT NULL,
    rangeType varchar(50) NOT NULL,
    rank varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableGenes;
CREATE TABLE actionableGenes
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    event varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableFusions;
CREATE TABLE actionableFusions
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    geneUp varchar(50) NOT NULL,
    minExonUp int,
    maxExonUp int,
    geneDown varchar(50) NOT NULL,
    minExonDown int,
    maxExonDown int,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableCharacteristics;
CREATE TABLE actionableCharacteristics
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    name varchar(50) NOT NULL,
    comparator varchar(50),
    cutOff int,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableHla;
CREATE TABLE actionableHla
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    HLAType varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(2000),
    treatment varchar(500) NOT NULL,
    drugClasses varchar(500),
    applicableCancerType varchar(100) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(500),
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS knownHotspots;
CREATE TABLE knownHotspots
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    chromosome varchar(50) NOT NULL,
    position varchar(50) NOT NULL,
    ref varchar(50) NOT NULL,
    alt varchar(50) NOT NULL,
    inputGene varchar(50) NOT NULL,
    inputTranscript varchar(50),
    inputProteinAnnotation varchar(50) NOT NULL,
    inputSource varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS knownCodons;
CREATE TABLE knownCodons
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    transcript varchar(50) NOT NULL,
    chromosome varchar(50) NOT NULL,
    start varchar(50) NOT NULL,
    end varchar(50) NOT NULL,
    mutationType varchar(50) NOT NULL,
    codonRank varchar(50) NOT NULL,
    sources varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS knownExons;
CREATE TABLE knownExons
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    transcript varchar(50) NOT NULL,
    chromosome varchar(50) NOT NULL,
    start varchar(50) NOT NULL,
    end varchar(50) NOT NULL,
    mutationType varchar(50) NOT NULL,
    exonRank varchar(50) NOT NULL,
    sources varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS knownCopyNumbers;
CREATE TABLE knownCopyNumbers
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(50) NOT NULL,
    type varchar(50) NOT NULL,
    sources varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS knownFusionPairs;
CREATE TABLE knownFusionPairs
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    geneUp varchar(50) NOT NULL,
    minExonUp varchar(50),
    maxExonUp varchar(50),
    geneDown varchar(50) NOT NULL,
    minExonDown varchar(50),
    maxExonDown varchar(50),
    sources varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actin;
CREATE TABLE actin
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    trial varchar(50) NOT NULL,
    cohort varchar(50),
    rule varchar(50) NOT NULL,
    gene varchar(50),
    mutation varchar(50),
    isUsedAsInclusion varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

SET FOREIGN_KEY_CHECKS = 1;