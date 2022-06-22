SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS actionableHotspots;
CREATE TABLE actionableHotspots
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    chromosome varchar(50) NOT NULL,
    position int NOT NULL,
    ref varchar(50) NOT NULL,
    alt varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
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
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
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
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableFusions;
CREATE TABLE actionableFusions
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    geneUp varchar(50) NOT NULL,
    minExonUp int NOT NULL,
    maxExonUp int NOT NULL,
    geneDown varchar(50) NOT NULL,
    minExonDown int NOT NULL,
    maxExonDown int NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableCharacteristics;
CREATE TABLE actionableCharacteristics
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    name varchar(50) NOT NULL,
    comparator varchar(50) NOT NULL,
    cutOff int NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actionableHla;
CREATE TABLE actionableHla
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    HLAType varchar(50) NOT NULL,
    source varchar(50) NOT NULL,
    sourceEvent varchar(50) NOT NULL,
    sourceUrls varchar(50) NOT NULL,
    treatment varchar(50) NOT NULL,
    drugClasses varchar(50) NOT NULL,
    applicableCancerType varchar(50) NOT NULL,
    applicableDoid varchar(50) NOT NULL,
    blacklistCancerTypes varchar(50) NOT NULL,
    level varchar(50) NOT NULL,
    direction varchar(50) NOT NULL,
    evidenceUrls varchar(50) NOT NULL,
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
    inputTranscript varchar(50) NOT NULL,
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
    minExonUp varchar(50) NOT NULL,
    maxExonUp varchar(50) NOT NULL,
    geneDown varchar(50) NOT NULL,
    minExonDown varchar(50) NOT NULL,
    maxExonDown varchar(50) NOT NULL,
    sources varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS iclusionStudy;
CREATE TABLE iclusionStudy
(   id int NOT NULL AUTO_INCREMENT,
    idDB varchar(50) NOT NULL,
    acronym varchar(50) NOT NULL,
    title varchar(50) NOT NULL,
    eudra varchar(50) NOT NULL,
    nct varchar(50) NOT NULL,
    ipn varchar(50) NOT NULL,
    ccmo varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS iclusionStudyTumorLocations;
CREATE TABLE iclusionStudyTumorLocations
(   id int NOT NULL AUTO_INCREMENT,
    tumorLocationId int NOT NULL,
    tumorLocation varchar(50) NOT NULL,
    FOREIGN KEY (tumorLocationId) REFERENCES iclusionStudy(id),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS iclusionStudyBlacklistedTumorLocations;
CREATE TABLE iclusionStudyBlacklistedTumorLocations
(   id int NOT NULL AUTO_INCREMENT,
    blacklistedTumorLocationId int NOT NULL,
    blacklistedTumorLocation varchar(50) NOT NULL,
    FOREIGN KEY (blacklistedTumorLocationId) REFERENCES iclusionStudy(id),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS iclusionStudyMutationConditions;
CREATE TABLE iclusionStudyMutationConditions
(   id int NOT NULL AUTO_INCREMENT,
    mutationConditionId int NOT NULL,
    gene varchar(50) NOT NULL,
    mutation varchar(50) NOT NULL,
    FOREIGN KEY (mutationConditionId) REFERENCES iclusionStudy(id),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS actin;
CREATE TABLE actin
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    trial varchar(50) NOT NULL,
    cohort varchar(50) NOT NULL,
    rule varchar(50) NOT NULL,
    gene varchar(50) NOT NULL,
    mutation varchar(50) NOT NULL,
    isUsedAsInclusion varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

SET FOREIGN_KEY_CHECKS = 1;