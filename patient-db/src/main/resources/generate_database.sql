SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS patient;
CREATE TABLE patient
(   id int NOT NULL AUTO_INCREMENT,
    cpctId varchar(255) DEFAULT NULL,
    registrationDate DATE,
    gender varchar(10),
    ethnicity varchar(255),
    hospital varchar(255),
    birthYear int,
    primaryTumorLocation varchar(255),
    deathDate DATE,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS sample;
CREATE TABLE sample
(   sampleId varchar(20) NOT NULL,
    patientId int NOT NULL,
    arrivalDate DATE NOT NULL,
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS biopsy;
CREATE TABLE biopsy
(   id int NOT NULL,
    sampleId varchar(20),
    patientId int NOT NULL,
    biopsyLocation varchar(255),
    biopsyDate DATE,
    PRIMARY KEY (id),
    FOREIGN KEY (sampleId) REFERENCES sample(sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS treatment;
CREATE TABLE treatment
 (  id int NOT NULL,
    biopsyId int,
    patientId int NOT NULL,
    treatmentGiven varchar(3),
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (biopsyId) REFERENCES biopsy(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
 );

SET FOREIGN_KEY_CHECKS = 1;

DROP TABLE IF EXISTS drug;
CREATE TABLE drug
 (  id int NOT NULL AUTO_INCREMENT,
    treatmentId int,
    patientId int NOT NULL,
    startDate DATE,
    endDate DATE,
    name varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentId) REFERENCES treatment(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
 );

DROP TABLE IF EXISTS treatmentResponse;
 CREATE TABLE treatmentResponse
  (  id int NOT NULL AUTO_INCREMENT,
     treatmentId int,
     patientId int NOT NULL,
     measurementDone varchar(5),
     responseDate DATE,
     response varchar(25),
     PRIMARY KEY (id),
     FOREIGN KEY (treatmentId) REFERENCES treatment(id),
     FOREIGN KEY (patientId) REFERENCES patient(id)
  );

DROP TABLE IF EXISTS somaticVariant;
CREATE TABLE somaticVariant
(   id int NOT NULL AUTO_INCREMENT,
    sampleId varchar(20) NOT NULL,
    patientId int NOT NULL,
    gene varchar(255) NOT NULL,
    position varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    cosmicId varchar(255),
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (sampleId) REFERENCES sample(sampleId),
    FOREIGN KEY (patientId) references patient(id)
);

DROP TABLE IF EXISTS comprehensiveSomaticVariant;
CREATE TABLE comprehensiveSomaticVariant
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    position int not null,
    filter varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    adjustedVaf DOUBLE PRECISION NOT NULL,
    adjustedCopyNumber DOUBLE PRECISION NOT NULL,
    highConfidence BOOLEAN NOT NULL,
    trinucleotideContext varchar(3) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(filter)
);

DROP TABLE IF EXISTS copyNumber;
CREATE TABLE copyNumber
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    ratioSupport BOOLEAN NOT NULL,
    structuralVariantSupport varchar(255) NOT NULL,
    bafCount int not null,
    observedBaf DOUBLE PRECISION not null,
    actualBaf DOUBLE PRECISION not null,
    copyNumber DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS purityRange;
CREATE TABLE purityRange
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    purity DOUBLE PRECISION not null,
    normFactor DOUBLE PRECISION not null,
    score DOUBLE PRECISION not null,
    ploidy DOUBLE PRECISION not null,
    diploidProportion DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS purityScore;
CREATE TABLE purityScore
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    polyclonalProportion DOUBLE PRECISION not null,
    minPurity DOUBLE PRECISION not null,
    maxPurity DOUBLE PRECISION not null,
    minPloidy DOUBLE PRECISION not null,
    maxPloidy DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

CREATE OR REPLACE VIEW purity AS
SELECT p.*, s.polyclonalProportion, s.minPurity, s.maxPurity, s.minPloidy, s.maxPloidy
FROM purityRange p, purityScore s
WHERE p.sampleId = s.sampleId
  AND (p.sampleId, p.score) IN (SELECT sampleId, MIN(score) FROM purityRange GROUP BY sampleId);

DROP TABLE IF EXISTS copyNumberRegion;
CREATE TABLE copyNumberRegion
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    ratioSupport BOOLEAN NOT NULL,
    structuralVariantSupport varchar(255) NOT NULL,
    bafCount int not null,
    observedBaf DOUBLE PRECISION not null,
    observedTumorRatio DOUBLE PRECISION not null,
    observedNormalRatio DOUBLE PRECISION not null,
    observedGCContent DOUBLE PRECISION not null,
    observedNonNPercentage DOUBLE PRECISION not null,
    observedMappablePercentage DOUBLE PRECISION not null,
    modelPloidy int not null,
    modelBaf DOUBLE PRECISION not null,
    modelTumorRatio DOUBLE PRECISION not null,
    actualTumorCopyNumber DOUBLE PRECISION not null,
    cnvDeviation DOUBLE PRECISION not null,
    bafDeviation DOUBLE PRECISION not null,
    highConfidenceBaf DOUBLE PRECISION not null,
    highConfidenceCopyNumber DOUBLE PRECISION not null,
    fittedBaf DOUBLE PRECISION not null,
    fittedCopyNumber DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS structuralVariant;
CREATE TABLE structuralVariant
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    startChromosome varchar(255) NOT NULL,
    endChromosome varchar(255) NOT NULL,
    startPosition int not null,
    endPosition int not null,
    startOrientation tinyint not null,
    endOrientation tinyint not null,
    startHomologySequence varchar(255) not null,
    endHomologySequence varchar(255) not null,
    insertSequence varchar(255) not null,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS geneCopyNumber;
CREATE TABLE geneCopyNumber
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(20) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    gene varchar(255) NOT NULL,
    chromosomeBand varchar(255) NOT NULL,
    transcriptId varchar(255) NOT NULL,
    transcriptVersion int not null,
    minCopyNumber DOUBLE PRECISION not null,
    maxCopyNumber DOUBLE PRECISION not null,
    meanCopyNumber DOUBLE PRECISION not null,
    regions int not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS clinicalLogs;
CREATE TABLE clinicalLogs
(   id int NOT NULL AUTO_INCREMENT,
    eventDate TIMESTAMP,
    level varchar(100),
    patientId varchar(20),
    ecrfItem varchar(100),
    formStatus varchar(5),
    formLocked varchar(5),
    message varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS ecrf;
CREATE TABLE ecrf
(   id int NOT NULL AUTO_INCREMENT,
    patientId varchar(20),
    studyEvent varchar(100),
    studyEventKey int not null,
    form varchar(100),
    formKey int not null,
    itemGroup varchar(100),
    itemGroupKey int not null,
    item varchar(100),
    itemValue varchar(1500),
    status varchar(5),
    locked varchar(5),
    sequenced varchar(5),
    PRIMARY KEY (id),
    INDEX(patientId),
    INDEX(studyEvent),
    INDEX(form),
    INDEX(itemGroup),
    INDEX(item),
    INDEX(itemValue)
);

DROP TABLE IF EXISTS formsMetadata;
CREATE TABLE formsMetadata
(   id int NOT NULL,
    tableName varchar(20),
    form varchar(20),
    status varchar(5),
    locked varchar(5),
    UNIQUE KEY (id, tableName, form)
);
