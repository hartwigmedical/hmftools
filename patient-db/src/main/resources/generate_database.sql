SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS patient;
CREATE TABLE patient
(   id int NOT NULL AUTO_INCREMENT,
    patientIdentifier varchar(50) UNIQUE,
    blacklisted BOOLEAN NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS baseline;
CREATE TABLE baseline
(   patientId int NOT NULL,
    registrationDate DATE,
    informedConsentDate DATE,
    gender varchar(10),
    hospital varchar(255),
    birthYear int,
    primaryTumorLocation varchar(255),
    primaryTumorSubLocation varchar(255),
    primaryTumorType varchar(255),
    primaryTumorSubType varchar(255),
    primaryTumorExtraDetails varchar(255),
    primaryTumorOverridden BOOLEAN,
    deathDate DATE,
    hasSystemicPreTreatment varchar(3),
    hasRadiotherapyPreTreatment varchar(3),
    preTreatments varchar(800),
    preTreatmentsType varchar(510),
    preTreatmentsMechanism varchar(510),
    PRIMARY KEY (patientId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS doidNode;
CREATE TABLE doidNode
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    doid varchar(255) NOT NULL,
    doidTerm varchar(255) NOT NULL,
    snomedConceptId varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS snomed;
CREATE TABLE snomed
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    snomedConceptId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS preTreatmentDrug;
CREATE TABLE preTreatmentDrug
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    startDate DATE,
    endDate DATE,
    name varchar(800),
    type varchar(255),
    mechanism varchar(255),
    bestResponse varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS sample;
CREATE TABLE sample
(   sampleId varchar(255) NOT NULL,
    patientId int NOT NULL,
    setName varchar(255) NOT NULL,
    arrivalDate DATE NOT NULL,
    samplingDate DATE,
    dnaNanograms int,
    limsPrimaryTumor varchar(255),
    pathologyTumorPercentage varchar(100),
    PRIMARY KEY (sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS rna;
CREATE TABLE rna
(   sampleId varchar(255) NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS snpcheck;
CREATE TABLE snpcheck
(   sampleId varchar(255) NOT NULL,
    isPass BOOLEAN NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS biopsy;
CREATE TABLE biopsy
(   id int NOT NULL,
    sampleId varchar(255),
    patientId int NOT NULL,
    biopsyTaken varchar(255),
    biopsyEvaluable varchar(255),
    biopsyType varchar(255),
    biopsySite varchar(255),
    biopsyLocation varchar(255),
    biopsyDate DATE,
    PRIMARY KEY (id),
    FOREIGN KEY (sampleId) REFERENCES sample(sampleId),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS treatment;
CREATE TABLE treatment
(   id int NOT NULL,
    biopsyId int,
    patientId int NOT NULL,
    treatmentGiven varchar(3),
    radiotherapyGiven varchar(3),
    startDate DATE,
    endDate DATE,
    name varchar(800),
    type varchar(255),
    mechanism varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (biopsyId) REFERENCES biopsy(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS drug;
CREATE TABLE drug
(   id int NOT NULL AUTO_INCREMENT,
    treatmentId int,
    patientId int NOT NULL,
    startDate DATE,
    endDate DATE,
    name varchar(800),
    type varchar(255),
    mechanism varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentId) REFERENCES treatment(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS treatmentResponse;
CREATE TABLE treatmentResponse
(   id int NOT NULL AUTO_INCREMENT,
    treatmentId int,
    patientId int NOT NULL,
    measurementDone varchar(5),
    boneOnlyDisease varchar(5),
    responseDate DATE,
    response varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (treatmentId) REFERENCES treatment(id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS tumorMarker;
CREATE TABLE tumorMarker
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    date DATE,
    marker varchar(50),
    measurement varchar(50),
    unit varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS ranoMeasurement;
CREATE TABLE ranoMeasurement
(   id int NOT NULL AUTO_INCREMENT,
    patientId int NOT NULL,
    responseDate DATE,
    therapyGiven varchar(50),
    targetLesionResponse varchar(50),
    noTargetLesionResponse varchar(50),
    overallResponse varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (patientId) REFERENCES patient(id)
);

DROP TABLE IF EXISTS clinicalFindings;
CREATE TABLE clinicalFindings
(   id int NOT NULL AUTO_INCREMENT,
    level varchar(30),
    patientId varchar(20),
    formStatus varchar(30),
    formLocked varchar(5),
    message varchar(1000),
    details varchar(1000),
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS cpctEcrf;
CREATE TABLE cpctEcrf
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
    status varchar(30),
    locked varchar(5),
    sequenced varchar(5),
    fieldName varchar(100),
    relevant varchar(5),
    PRIMARY KEY (id),
    INDEX(patientId),
    INDEX(studyEvent),
    INDEX(form),
    INDEX(itemGroup),
    INDEX(item),
    INDEX(itemValue (255)),
    INDEX(status),
    INDEX(locked),
    INDEX(sequenced),
    INDEX(fieldName),
    INDEX(relevant)
);

DROP TABLE IF EXISTS cpctEcrfDatamodel;
CREATE TABLE cpctEcrfDatamodel
(   fieldName varchar(100),
    description varchar(500),
    codeList varchar(3000),
    relevant varchar(5)
);

DROP TABLE IF EXISTS drupEcrf;
CREATE TABLE drupEcrf
(   id int NOT NULL AUTO_INCREMENT,
    patientId varchar(20),
    studyEvent varchar(100),
    studyEventKey int not null,
    form varchar(100),
    formKey int not null,
    itemGroup varchar(100),
    itemGroupKey int not null,
    item varchar(100),
    itemValue varchar(15000),
    status varchar(30),
    locked varchar(5),
    sequenced varchar(5),
    fieldName varchar(100),
    relevant varchar(5),
    PRIMARY KEY (id),
    INDEX(patientId),
    INDEX(studyEvent),
    INDEX(form),
    INDEX(itemGroup),
    INDEX(item),
    INDEX(status),
    INDEX(sequenced),
    INDEX(fieldName),
    INDEX(relevant)
);

DROP TABLE IF EXISTS drupEcrfDatamodel;
CREATE TABLE drupEcrfDatamodel
(   fieldName varchar(100),
    description varchar(500),
    codeList varchar(5000),
    relevant varchar(5)
);

DROP TABLE IF EXISTS formsMetadata;
CREATE TABLE formsMetadata
(   id int NOT NULL,
    tableName varchar(20),
    form varchar(20),
    status varchar(30),
    locked varchar(5),
    UNIQUE KEY (id, tableName, form)
);

DROP TABLE if EXISTS metric;
CREATE TABLE metric
(   sampleId varchar(255) NOT NULL,
    refMeanCoverage DOUBLE PRECISION NOT NULL,
    refSdCoverage DOUBLE PRECISION NOT NULL,
    refMedianCoverage int NOT NULL,
    refMadCoverage int NOT NULL,
    refPctExcAdapter DOUBLE PRECISION,
    refPctExcMapQ DOUBLE PRECISION NOT NULL,
    refPctExcDupe DOUBLE PRECISION NOT NULL,
    refPctExcUnpaired DOUBLE PRECISION NOT NULL,
    refPctExcBaseQ DOUBLE PRECISION NOT NULL,
    refPctExcOverlap DOUBLE PRECISION NOT NULL,
    refPctExcCapped DOUBLE PRECISION NOT NULL,
    refPctExcTotal DOUBLE PRECISION NOT NULL,
    refCoverage1xPercentage DOUBLE PRECISION,
    refCoverage10xPercentage DOUBLE PRECISION NOT NULL,
    refCoverage20xPercentage DOUBLE PRECISION NOT NULL,
    refCoverage30xPercentage DOUBLE PRECISION NOT NULL,
    refCoverage60xPercentage DOUBLE PRECISION NOT NULL,
    tumorMeanCoverage DOUBLE PRECISION NOT NULL,
    tumorSdCoverage DOUBLE PRECISION NOT NULL,
    tumorMedianCoverage int NOT NULL,
    tumorMadCoverage int NOT NULL,
    tumorPctExcAdapter DOUBLE PRECISION,
    tumorPctExcMapQ DOUBLE PRECISION NOT NULL,
    tumorPctExcDupe DOUBLE PRECISION NOT NULL,
    tumorPctExcUnpaired DOUBLE PRECISION NOT NULL,
    tumorPctExcBaseQ DOUBLE PRECISION NOT NULL,
    tumorPctExcOverlap DOUBLE PRECISION NOT NULL,
    tumorPctExcCapped DOUBLE PRECISION NOT NULL,
    tumorPctExcTotal DOUBLE PRECISION NOT NULL,
    tumorCoverage1xPercentage DOUBLE PRECISION,
    tumorCoverage10xPercentage DOUBLE PRECISION NOT NULL,
    tumorCoverage20xPercentage DOUBLE PRECISION NOT NULL,
    tumorCoverage30xPercentage DOUBLE PRECISION NOT NULL,
    tumorCoverage60xPercentage DOUBLE PRECISION NOT NULL,
    sufficientCoverage BOOLEAN NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE if EXISTS flagstat;
CREATE TABLE flagstat
(   sampleId varchar(255) NOT NULL,
    refUniqueReadCount BIGINT NOT NULL,
    refSecondaryCount BIGINT NOT NULL,
    refSupplementaryCount BIGINT NOT NULL,
    refDuplicateProportion DOUBLE PRECISION NOT NULL,
    refMappedProportion DOUBLE PRECISION NOT NULL,
    refPairedInSequencingProportion DOUBLE PRECISION NOT NULL,
    refProperlyPairedProportion DOUBLE PRECISION NOT NULL,
    refWithItselfAndMateMappedProportion DOUBLE PRECISION NOT NULL,
    refSingletonProportion DOUBLE PRECISION NOT NULL,
    tumorUniqueReadCount BIGINT NOT NULL,
    tumorSecondaryCount BIGINT NOT NULL,
    tumorSupplementaryCount BIGINT NOT NULL,
    tumorDuplicateProportion DOUBLE PRECISION NOT NULL,
    tumorMappedProportion DOUBLE PRECISION NOT NULL,
    tumorPairedInSequencingProportion DOUBLE PRECISION NOT NULL,
    tumorProperlyPairedProportion DOUBLE PRECISION NOT NULL,
    tumorWithItselfAndMateMappedProportion DOUBLE PRECISION NOT NULL,
    tumorSingletonProportion DOUBLE PRECISION NOT NULL,
    passQC BOOLEAN NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS somaticVariant;
CREATE TABLE somaticVariant
(   id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    position int not null,
    filter varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    genesEffected int not null,
    worstEffectTranscript varchar(255) NOT NULL,
    worstEffect varchar(255) NOT NULL,
    worstCodingEffect varchar(255) NOT NULL,
    canonicalEffect varchar(255) NOT NULL,
    canonicalCodingEffect varchar(255) NOT NULL,
    canonicalHgvsCodingImpact varchar(255) NOT NULL,
    canonicalHgvsProteinImpact varchar(255) NOT NULL,
    microhomology varchar(255) NOT NULL,
    repeatSequence varchar(255) NOT NULL,
    repeatCount int NOT NULL,
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    adjustedVaf DOUBLE PRECISION NOT NULL,
    variantCopyNumber DOUBLE PRECISION NOT NULL,
    copyNumber DOUBLE PRECISION NOT NULL,
    tier varchar(20) NOT NULL,
    trinucleotideContext varchar(3) NOT NULL,
    subclonalLikelihood DOUBLE PRECISION NOT NULL,
    biallelic BOOLEAN NOT NULL,
    hotspot varchar(20) NOT NULL,
    mappability DOUBLE PRECISION NOT NULL,
    germlineStatus varchar(255) NOT NULL,
    minorAlleleCopyNumber DOUBLE PRECISION NOT NULL,
    recovered BOOLEAN NOT NULL,
    kataegis varchar(20) NOT NULL,
    referenceAlleleReadCount int,
    referenceTotalReadCount int,
    rnaAlleleReadCount int,
    rnaTotalReadCount int,
    localPhaseSet int,
    localRealignmentSet int,
    phasedInframeIndel int,
    qual double precision not null,
    reported BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(filter),
    INDEX(type),
    INDEX(gene)
);

DROP TABLE IF EXISTS purity;
CREATE TABLE purity
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    version VARCHAR(255) NOT NULL,
    sampleId varchar(255) NOT NULL,
    gender varchar(255) NOT NULL,
    fitMethod varchar(255) NOT NULL,
    qcStatus varchar(255) NOT NULL,
    purity DOUBLE PRECISION not null,
    normFactor DOUBLE PRECISION not null,
    score DOUBLE PRECISION not null,
    somaticPenalty DOUBLE PRECISION not null,
    ploidy DOUBLE PRECISION not null,
    diploidProportion DOUBLE PRECISION not null,
    polyclonalProportion DOUBLE PRECISION not null,
    wholeGenomeDuplication BOOLEAN NOT NULL,
    minPurity DOUBLE PRECISION not null,
    maxPurity DOUBLE PRECISION not null,
    minPloidy DOUBLE PRECISION not null,
    maxPloidy DOUBLE PRECISION not null,
    minDiploidProportion DOUBLE PRECISION not null,
    maxDiploidProportion DOUBLE PRECISION not null,
    msIndelsPerMb DOUBLE PRECISION not null,
    msStatus varchar(10) not null,
    tmbPerMb DOUBLE PRECISION not null,
    tmbStatus varchar(10) not null,
    tml INT not null,
    tmlStatus varchar(10) not null,
    svTmb INT not null DEFAULT 0,
    deletedGenes INT not null DEFAULT 0,
    copyNumberSegments INT not null DEFAULT 0,
    unsupportedCopyNumberSegments INT not null DEFAULT 0,
    contamination DOUBLE PRECISION not null DEFAULT 0,
    germlineAberration varchar(255) not null DEFAULT "NONE",
    amberGender varchar(255) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS purityRange;
CREATE TABLE purityRange
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    purity DOUBLE PRECISION not null,
    normFactor DOUBLE PRECISION not null,
    score DOUBLE PRECISION not null,
    somaticPenalty DOUBLE PRECISION not null,
    ploidy DOUBLE PRECISION not null,
    diploidProportion DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS copyNumber;
CREATE TABLE copyNumber
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    segmentStartSupport varchar(255) NOT NULL,
    segmentEndSupport varchar(255) NOT NULL,
    depthWindowCount int not null,
    bafCount int not null,
    observedBaf DOUBLE PRECISION not null,
    baf DOUBLE PRECISION not null,
    copyNumber DOUBLE PRECISION not null,
    minorAlleleCopyNumber DOUBLE PRECISION not null,
    majorAlleleCopyNumber DOUBLE PRECISION not null,
    copyNumberMethod varchar(255) NOT NULL,
    gcContent DOUBLE PRECISION not null,
    minStart int not null,
    maxStart int not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS copyNumberGermline;
CREATE TABLE copyNumberGermline
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    segmentStartSupport varchar(255) NOT NULL,
    segmentEndSupport varchar(255) NOT NULL,
    depthWindowCount int not null,
    bafCount int not null,
    observedBaf DOUBLE PRECISION not null,
    baf DOUBLE PRECISION not null,
    copyNumber DOUBLE PRECISION not null,
    minorAlleleCopyNumber DOUBLE PRECISION not null,
    majorAlleleCopyNumber DOUBLE PRECISION not null,
    copyNumberMethod varchar(255) NOT NULL,
    gcContent DOUBLE PRECISION not null,
    minStart int not null,
    maxStart int not null,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS geneCopyNumber;
CREATE TABLE geneCopyNumber
(   id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    start int not null,
    end int not null,
    gene varchar(255) NOT NULL,
    chromosomeBand varchar(255) NOT NULL,
    transcriptId varchar(255) NOT NULL,
    transcriptVersion int not null,
    minCopyNumber DOUBLE PRECISION not null,
    maxCopyNumber DOUBLE PRECISION not null,
    somaticRegions int not null,
    germlineHomDeletionRegions int not null,
    germlineHetToHomDeletionRegions int not null,
    minRegions int not null,
    minRegionStart int not null,
    minRegionEnd int not null,
    minRegionStartSupport varchar(255) NOT NULL,
    minRegionEndSupport varchar(255) NOT NULL,
    minRegionMethod varchar(255) NOT NULL,
    minMinorAlleleCopyNumber DOUBLE PRECISION not null,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);

DROP TABLE IF EXISTS structuralVariant;
CREATE TABLE structuralVariant
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    svId INT NOT NULL,
    sampleId varchar(255) NOT NULL,
    startChromosome varchar(255) NOT NULL,
    endChromosome varchar(255),
    startPosition int not null,
    endPosition int,
    startOrientation tinyint not null,
    endOrientation tinyint,
    startHomologySequence varchar(255) not null,
    endHomologySequence varchar(255),
    startAF DOUBLE PRECISION,
    endAF DOUBLE PRECISION,
    junctionCopyNumber DOUBLE PRECISION,
    adjustedAFStart DOUBLE PRECISION,
    adjustedAFEnd DOUBLE PRECISION,
    adjustedCopyNumberStart DOUBLE PRECISION,
    adjustedCopyNumberEnd DOUBLE PRECISION,
    adjustedCopyNumberChangeStart DOUBLE PRECISION,
    adjustedCopyNumberChangeEnd DOUBLE PRECISION,
    insertSequence varchar(2048) not null,
    type varchar(255) NOT NULL,
    filter varchar(255) NOT NULL,
    imprecise BOOLEAN NOT NULL,
    qualScore DOUBLE PRECISION,
    event varchar(255),
    startTumorVariantFragmentCount int,
    startTumorReferenceFragmentCount int,
    startNormalVariantFragmentCount int,
    startNormalReferenceFragmentCount int,
    endTumorVariantFragmentCount int,
    endTumorReferenceFragmentCount int,
    endNormalVariantFragmentCount int,
    endNormalReferenceFragmentCount int,
    startIntervalOffsetStart int,
    startIntervalOffsetEnd int,
    endIntervalOffsetStart int,
    endIntervalOffsetEnd int,
    inexactHomologyOffsetStart int,
    inexactHomologyOffsetEnd int,
    startLinkedBy varchar(1024),
    endLinkedBy varchar(1024),
    vcfId varchar(255),
    recovered BOOLEAN NOT NULL,
    recoveryMethod varchar(64),
    recoveryFilter varchar(255),
    startRefContext varchar(255),
    endRefContext varchar(255),
    insertSequenceAlignments varchar(512),
    insertSequenceRepeatClass varchar(64),
    insertSequenceRepeatType varchar(64),
    insertSequenceRepeatOrientation tinyint,
    insertSequenceRepeatCoverage DOUBLE PRECISION,
    startAnchoringSupportDistance int,
    endAnchoringSupportDistance int,
    PRIMARY KEY (id),
    INDEX(sampleId, svId)
);

DROP TABLE IF EXISTS svAnnotation;
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
    junctionCopyNumberMin DOUBLE PRECISION,
    junctionCopyNumberMax DOUBLE PRECISION,
    geneStart VARCHAR(20),
    geneEnd varchar(20),
    replicationTimingStart DOUBLE PRECISION,
    replicationTimingEnd DOUBLE PRECISION,
    localTopologyIdStart INT,
    localTopologyIdEnd INT,
    localTopologyStart varchar(20),
    localTopologyEnd VARCHAR(20),
    localTICountStart INT,
    localTICountEnd INT,
    PRIMARY KEY (id),
    INDEX(sampleId, svId),
    INDEX(clusterId),
    INDEX(svId)
);

DROP TABLE IF EXISTS svCluster;
CREATE TABLE svCluster
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NOT NULL,
    category VARCHAR(20),
    synthetic BOOLEAN NOT NULL,
    resolvedType VARCHAR(20),
    clusterCount INT,
    clusterDesc VARCHAR(50),
    PRIMARY KEY (id),
    INDEX(sampleId, clusterId),
    INDEX(clusterId)
);

DROP TABLE IF EXISTS svLink;
CREATE TABLE svLink
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NOT NULL,
    chainId INT NOT NULL,
    chainIndex VARCHAR(50) NOT NULL,
    chainLinkCount INT NOT NULL,
    lowerSvId INT NOT NULL,
    upperSvId INT NOT NULL,
    lowerBreakendIsStart BOOLEAN NOT NULL,
    upperBreakendIsStart BOOLEAN NOT NULL,
    chromosome VARCHAR(10),
    arm VARCHAR(3),
    assembled BOOLEAN NOT NULL,
    traversedSVCount INT,
    linkLength INT,
    junctionCopyNumber DOUBLE PRECISION,
    junctionCopyNumberUncertainty DOUBLE PRECISION,
    pseudogeneInfo varchar(255),
    ecDna BOOLEAN,
    PRIMARY KEY (id),
    INDEX(sampleId, clusterId),
    INDEX(clusterId)
);

DROP TABLE IF EXISTS svDriver;
CREATE TABLE svDriver
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    clusterId INT NULL,
    gene VARCHAR(50) NOT NULL,
    eventType VARCHAR(50),
    PRIMARY KEY (id),
    INDEX(sampleId, clusterId),
    INDEX(clusterId)
);

DROP TABLE IF EXISTS svBreakend;
CREATE TABLE svBreakend
(   id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
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
    codingContext VARCHAR(20),
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

DROP TABLE IF EXISTS svFusion;
CREATE TABLE svFusion
(   id INT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    fivePrimeBreakendId INT UNSIGNED NOT NULL,
    threePrimeBreakendId INT UNSIGNED NOT NULL,
    name VARCHAR(50) NOT NULL,
    reported BOOLEAN NOT NULL,
    reportedType varchar(50) NOT NULL,
    phased VARCHAR(20) NOT NULL,
    likelihood VARCHAR(10) NOT NULL,
    chainLength INT,
    chainLinks INT,
    chainTerminated BOOLEAN,
    domainsKept VARCHAR(255),
    domainsLost VARCHAR(255),
    skippedExonsUp INT,
    skippedExonsDown INT,
    fusedExonUp INT,
    fusedExonDown INT,
    PRIMARY KEY (id),
    INDEX(fivePrimeBreakendId),
    INDEX(threePrimeBreakendId),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS viralInsertion;
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

DROP TABLE IF EXISTS signature;
CREATE TABLE signature
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    signature VARCHAR(50) NOT NULL,
    allocation DOUBLE PRECISION,
    percent DOUBLE PRECISION,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS canonicalTranscript;
CREATE TABLE canonicalTranscript
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    assembly varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    geneId varchar(255) NOT NULL,
    chromosomeBand varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    geneStart int not null,
    geneEnd int not null,
    transcriptId varchar(255) NOT NULL,
    transcriptVersion int not null,
    transcriptStart int not null,
    transcriptEnd int not null,
    exons int not null,
    exonStart int not null,
    exonEnd int not null,
    exonBases int not null,
    strand varchar(255) not null,
    codingStart int not null,
    codingEnd int not null,
    codingBases int not null,
    PRIMARY KEY (id),
    INDEX(gene),
    INDEX(transcriptId)
);

DROP TABLE IF EXISTS germlineVariant;
CREATE TABLE germlineVariant
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    position int not null,
    filter varchar(255) NOT NULL,
    refStatus varchar(20) NOT NULL,
    reported BOOLEAN NOT NULL,
    pathogenic varchar(50) NOT NULL,
    clinvarInfo VARCHAR(255),
    type varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    effect varchar(255) NOT NULL,
    codingEffect varchar(255) NOT NULL,
    hgvsCoding varchar(255) NOT NULL,
    hgvsProtein varchar(255) NOT NULL,
    microhomology varchar(255) NOT NULL,
    repeatSequence varchar(255) NOT NULL,
    repeatCount int NOT NULL,
    trinucleotideContext varchar(3) NOT NULL,
    alleleReadCount int NOT NULL,
    totalReadCount int NOT NULL,
    adjustedCopyNumber DOUBLE PRECISION NOT NULL,
    minorAlleleCopyNumber DOUBLE PRECISION NOT NULL,
    adjustedVaf DOUBLE PRECISION NOT NULL,
    biallelic BOOLEAN NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS germlineVariant2;
CREATE TABLE germlineVariant2
(   id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    position int not null,
    filter varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,

    ## SAGE
    qual double precision not null,
    tier varchar(20) NOT NULL,
    germlineGenotype varchar(255) NOT NULL,
    germlineAlleleReadCount int,
    germlineTotalReadCount int,
    rnaAlleleReadCount int,
    rnaTotalReadCount int,
    tumorAlleleReadCount int NOT NULL,
    tumorTotalReadCount int NOT NULL,
    localPhaseSet int,

    ### PURPLE ENRICHMENT
    adjustedVaf DOUBLE PRECISION NOT NULL,
    variantCopyNumber DOUBLE PRECISION NOT NULL,
    copyNumber DOUBLE PRECISION NOT NULL,
    biallelic BOOLEAN NOT NULL,
    minorAlleleCopyNumber DOUBLE PRECISION NOT NULL,

    ### PATHOGENIC
    clinvarInfo varchar(255) NOT NULL,
    pathogenicity varchar(255) NOT NULL,
    pathogenic BOOLEAN NOT NULL,

    ## SNP EFF ENRICHMENT
    gene varchar(255) NOT NULL,
    genesEffected int not null,
    worstEffectTranscript varchar(255) NOT NULL,
    worstEffect varchar(255) NOT NULL,
    worstCodingEffect varchar(255) NOT NULL,
    canonicalEffect varchar(255) NOT NULL,
    canonicalCodingEffect varchar(255) NOT NULL,
    canonicalHgvsCodingImpact varchar(255) NOT NULL,
    canonicalHgvsProteinImpact varchar(255) NOT NULL,

    ### REF GENOME ENRICHMENT
    microhomology varchar(255) NOT NULL,
    repeatSequence varchar(255) NOT NULL,
    repeatCount int NOT NULL,
    trinucleotideContext varchar(3) NOT NULL,

    ### DRIVER CATALOG
    hotspot varchar(20) NOT NULL,
    mappability DOUBLE PRECISION NOT NULL,
    reported BOOLEAN NOT NULL,

    ### NOT USED
    #recovered BOOLEAN NOT NULL,
    #germlineStatus varchar(255) NOT NULL,
    #kataegis varchar(20) NOT NULL,
    #localRealignmentSet int,
    #phasedInframeIndel int,
    #subclonalLikelihood DOUBLE PRECISION NOT NULL,

    PRIMARY KEY (id),
    INDEX(sampleId),
    INDEX(filter),
    INDEX(type),
    INDEX(gene)
);

DROP TABLE IF EXISTS driverCatalog;
CREATE TABLE driverCatalog
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    chromosomeBand varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    category varchar(255) NOT NULL,
    driver varchar(255) NOT NULL,
    likelihoodMethod varchar(255) NOT NULL,
    driverLikelihood DOUBLE PRECISION NOT NULL,
    missense int NOT NULL,
    nonsense int NOT NULL,
    splice int NOT NULL,
    frameshift int NOT NULL,
    inframe int NOT NULL,
    biallelic BOOLEAN NOT NULL,
    minCopyNumber DOUBLE PRECISION NOT NULL,
    maxCopyNumber DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);

DROP TABLE IF EXISTS chord;
CREATE TABLE chord
(   sampleId varchar(255) NOT NULL,
    BRCA1 DOUBLE PRECISION NOT NULL,
    BRCA2 DOUBLE PRECISION NOT NULL,
    hrd DOUBLE PRECISION NOT NULL,
    hrStatus varchar(255) NOT NULL,
    hrdType varchar(255) NOT NULL,
    remarksHrStatus varchar(255) NOT NULL,
    remarksHrdType varchar(255) NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS protect;
CREATE TABLE protect
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    event varchar(255) NOT NULL,
    germline BOOLEAN NOT NULL,
    reported BOOLEAN NOT NULL,
    treatment varchar(255) NOT NULL,
    onLabel BOOLEAN NOT NULL,
    level varchar(255) NOT NULL,
    direction varchar(255) NOT NULL,
    sources varchar(255) NOT NULL,
    urls varchar(2500) NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS peachCalls;
CREATE TABLE peachCalls
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    positionV37 varchar(255) NOT NULL,
    positionV38 varchar(255) NOT NULL,
    refV37 varchar(255) NOT NULL,
    refV38 varchar(255) NOT NULL,
    allele1 varchar(255) NOT NULL,
    allele2 varchar(255) NOT NULL,
    rsid varchar(255) NOT NULL,
    variantAnnotationV37 varchar(255) NOT NULL,
    filterV37 varchar(255) NOT NULL,
    variantAnnotationV38 varchar(255) NOT NULL,
    filterV38 varchar(255) NOT NULL,
    panelVersion varchar(255) NOT NULL,
    repoVersion varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS virusBreakend;
CREATE TABLE virusBreakend
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    taxidGenus int NOT NULL,
    nameGenus varchar(255) NOT NULL,
    readsGenusTree int NOT NULL,
    taxidSpecies int NOT NULL,
    nameSpecies varchar(255) NOT NULL,
    readsSpeciesTree int NOT NULL,
    taxidAssigned int NOT NULL,
    nameAssigned varchar(255) NOT NULL,
    readsAssignedTree int NOT NULL,
    readsAssignedDirect int NOT NULL,
    reference varchar(255) NOT NULL,
    referenceTaxid int NOT NULL,
    referenceKmerCount int NOT NULL,
    alternateKmerCount int NOT NULL,
    RName varchar(255) NOT NULL,
    startPos int NOT NULL,
    endPos int NOT NULL,
    numReads int NOT NULL,
    covBases int NOT NULL,
    coverage int NOT NULL,
    meanDepth int NOT NULL,
    meanBaseQ int NOT NULL,
    meanMapQ int NOT NULL,
    integrations int NOT NULL,
    qcStatus varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS peachGenotype;
CREATE TABLE peachGenotype
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    haplotype varchar(255) NOT NULL,
    `function` varchar(255) NOT NULL,
    linkedDrugs varchar(255) NOT NULL,
    urlPrescriptionInfo varchar(255) NOT NULL,
    panelVersion varchar(255) NOT NULL,
    repoVersion varchar(255) NOT NULL,
    PRIMARY KEY (id)
);

DROP TABLE IF EXISTS amberPatient;
CREATE TABLE amberPatient
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    sampleId varchar(255) NOT NULL,
    patientId int NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS amberAnonymous;
CREATE TABLE amberAnonymous
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    sampleId varchar(255) NOT NULL,
    hmfSampleId varchar(255) NOT NULL,
    deleted BOOLEAN NOT NULL,
    PRIMARY KEY (sampleId),
    INDEX(hmfSampleId)
);

DROP TABLE IF EXISTS amberMapping;
CREATE TABLE amberMapping
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    firstSampleId varchar(255) NOT NULL,
    secondSampleId varchar(255) NOT NULL,
    matches int NOT NULL,
    sites int NOT NULL,
    likelihood DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (firstSampleId, secondSampleId)
);

DROP TABLE IF EXISTS amberSample;
CREATE TABLE amberSample
(
    modified DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    sampleId varchar(255) NOT NULL,
    site1 TINYINT NOT NULL,
    site2 TINYINT NOT NULL,
    site3 TINYINT NOT NULL,
    site4 TINYINT NOT NULL,
    site5 TINYINT NOT NULL,
    site6 TINYINT NOT NULL,
    site7 TINYINT NOT NULL,
    site8 TINYINT NOT NULL,
    site9 TINYINT NOT NULL,
    site10 TINYINT NOT NULL,
    site11 TINYINT NOT NULL,
    site12 TINYINT NOT NULL,
    site13 TINYINT NOT NULL,
    site14 TINYINT NOT NULL,
    site15 TINYINT NOT NULL,
    site16 TINYINT NOT NULL,
    site17 TINYINT NOT NULL,
    site18 TINYINT NOT NULL,
    site19 TINYINT NOT NULL,
    site20 TINYINT NOT NULL,
    site21 TINYINT NOT NULL,
    site22 TINYINT NOT NULL,
    site23 TINYINT NOT NULL,
    site24 TINYINT NOT NULL,
    site25 TINYINT NOT NULL,
    site26 TINYINT NOT NULL,
    site27 TINYINT NOT NULL,
    site28 TINYINT NOT NULL,
    site29 TINYINT NOT NULL,
    site30 TINYINT NOT NULL,
    site31 TINYINT NOT NULL,
    site32 TINYINT NOT NULL,
    site33 TINYINT NOT NULL,
    site34 TINYINT NOT NULL,
    site35 TINYINT NOT NULL,
    site36 TINYINT NOT NULL,
    site37 TINYINT NOT NULL,
    site38 TINYINT NOT NULL,
    site39 TINYINT NOT NULL,
    site40 TINYINT NOT NULL,
    site41 TINYINT NOT NULL,
    site42 TINYINT NOT NULL,
    site43 TINYINT NOT NULL,
    site44 TINYINT NOT NULL,
    site45 TINYINT NOT NULL,
    site46 TINYINT NOT NULL,
    site47 TINYINT NOT NULL,
    site48 TINYINT NOT NULL,
    site49 TINYINT NOT NULL,
    site50 TINYINT NOT NULL,
    site51 TINYINT NOT NULL,
    site52 TINYINT NOT NULL,
    site53 TINYINT NOT NULL,
    site54 TINYINT NOT NULL,
    site55 TINYINT NOT NULL,
    site56 TINYINT NOT NULL,
    site57 TINYINT NOT NULL,
    site58 TINYINT NOT NULL,
    site59 TINYINT NOT NULL,
    site60 TINYINT NOT NULL,
    site61 TINYINT NOT NULL,
    site62 TINYINT NOT NULL,
    site63 TINYINT NOT NULL,
    site64 TINYINT NOT NULL,
    site65 TINYINT NOT NULL,
    site66 TINYINT NOT NULL,
    site67 TINYINT NOT NULL,
    site68 TINYINT NOT NULL,
    site69 TINYINT NOT NULL,
    site70 TINYINT NOT NULL,
    site71 TINYINT NOT NULL,
    site72 TINYINT NOT NULL,
    site73 TINYINT NOT NULL,
    site74 TINYINT NOT NULL,
    site75 TINYINT NOT NULL,
    site76 TINYINT NOT NULL,
    site77 TINYINT NOT NULL,
    site78 TINYINT NOT NULL,
    site79 TINYINT NOT NULL,
    site80 TINYINT NOT NULL,
    site81 TINYINT NOT NULL,
    site82 TINYINT NOT NULL,
    site83 TINYINT NOT NULL,
    site84 TINYINT NOT NULL,
    site85 TINYINT NOT NULL,
    site86 TINYINT NOT NULL,
    site87 TINYINT NOT NULL,
    site88 TINYINT NOT NULL,
    site89 TINYINT NOT NULL,
    site90 TINYINT NOT NULL,
    site91 TINYINT NOT NULL,
    site92 TINYINT NOT NULL,
    site93 TINYINT NOT NULL,
    site94 TINYINT NOT NULL,
    site95 TINYINT NOT NULL,
    site96 TINYINT NOT NULL,
    site97 TINYINT NOT NULL,
    site98 TINYINT NOT NULL,
    site99 TINYINT NOT NULL,
    site100 TINYINT NOT NULL,
    PRIMARY KEY (sampleId)
);

DROP TABLE IF EXISTS driverGenePanel;
CREATE TABLE driverGenePanel
(
    modified DATETIME NOT NULL,
    gene varchar(255) NOT NULL,
    reportMissense BOOLEAN NOT NULL,
    reportNonsense BOOLEAN NOT NULL,
    reportSplice BOOLEAN NOT NULL,
    reportDeletion BOOLEAN NOT NULL,
    reportDisruption BOOLEAN NOT NULL,
    reportAmplification BOOLEAN NOT NULL,
    reportSomaticHotspot BOOLEAN NOT NULL,
    reportGermlineVariant BOOLEAN NOT NULL,
    reportGermlineHotspot BOOLEAN NOT NULL,
    likelihoodType varchar(255) NOT NULL,
    PRIMARY KEY (gene)
);

DROP TABLE IF EXISTS rnaStatistics;
CREATE TABLE rnaStatistics
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    readLength int NOT NULL,
    totalFragments int NOT NULL,
    duplicates int NOT NULL,
    splicedPercent DOUBLE PRECISION NOT NULL,
    unsplicedPercent DOUBLE PRECISION NOT NULL,
    alternateSplicePercent DOUBLE PRECISION NOT NULL,
    chimericPercent DOUBLE PRECISION NOT NULL,
    fragmentLengthPct05 DOUBLE PRECISION NOT NULL,
    fragmentLengthPct50 DOUBLE PRECISION NOT NULL,
    fragmentLengthPct95 DOUBLE PRECISION NOT NULL,
    enrichedGenePercent DOUBLE PRECISION NOT NULL,
    medianGCRatio DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId)
);

DROP TABLE IF EXISTS geneExpression;
CREATE TABLE geneExpression
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene VARCHAR(20) NOT NULL,
    tpm DOUBLE PRECISION NOT NULL,
    splicedFragments int NOT NULL,
    unsplicedFragments int NOT NULL,
    medianTpmCancer DOUBLE PRECISION NOT NULL,
    percentileCancer DOUBLE PRECISION NOT NULL,
    medianTpmCohort DOUBLE PRECISION NOT NULL,
    percentileCohort DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);

DROP TABLE IF EXISTS novelSpliceJunction;
CREATE TABLE novelSpliceJunction
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    gene VARCHAR(20) NOT NULL,
    chromosome VARCHAR(10) NOT NULL,
    junctionStart int NOT NULL,
    junctionEnd int NOT NULL,
    type VARCHAR(20) NOT NULL,
    fragmentCount int NOT NULL,
    depthStart int NOT NULL,
    depthEnd int NOT NULL,
    regionStart VARCHAR(20) NOT NULL,
    regionEnd VARCHAR(20) NOT NULL,
    basesStart VARCHAR(20) NOT NULL,
    basesEnd VARCHAR(20) NOT NULL,
    cohortFrequency int NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, gene),
    INDEX(gene)
);

DROP TABLE IF EXISTS rnaFusion;
CREATE TABLE rnaFusion
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    name VARCHAR(50) NOT NULL,
    chromosomeUp VARCHAR(10) NOT NULL,
    chromosomeDown VARCHAR(10) NOT NULL,
    positionUp int NOT NULL,
    positionDown int NOT NULL,
    orientationUp tinyint NOT NULL,
    orientationDown tinyint NOT NULL,
    junctionTypeUp VARCHAR(20) NOT NULL,
    junctionTypeDown VARCHAR(20) NOT NULL,
    svType VARCHAR(5) NOT NULL,
    splitFragments int NOT NULL,
    realignedFrags int NOT NULL,
    discordantFrags int NOT NULL,
    depthUp int NOT NULL,
    depthDown int NOT NULL,
    maxAnchorLengthUp int NOT NULL,
    maxAnchorLengthDown int NOT NULL,
    cohortFrequency int NOT NULL,
    PRIMARY KEY (id),
    INDEX(sampleId, name),
    INDEX(name)
);

DROP TABLE IF EXISTS hlaType;
CREATE TABLE hlaType
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    status varchar(255) NOT NULL,
    typeA1 varchar(255) NOT NULL,
    typeA2 varchar(255) NOT NULL,
    typeB1 varchar(255) NOT NULL,
    typeB2 varchar(255) NOT NULL,
    typeC1 varchar(255) NOT NULL,
    typeC2 varchar(255) NOT NULL,
    somaticVariants int NOT NULL,
    PRIMARY KEY (id),
    KEY(sampleId)
);

DROP TABLE IF EXISTS hlaTypeDetails;
CREATE TABLE hlaTypeDetails
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    refUniqueCoverage INT NOT NULL,
    refSharedCoverage INT NOT NULL,
    refWildcardCoverage INT NOT NULL,
    tumorUniqueCoverage INT NOT NULL,
    tumorSharedCoverage INT NOT NULL,
    tumorWildcardCoverage INT NOT NULL,
    tumorCopyNumber DOUBLE PRECISION NOT NULL,
    somaticInframeIndel DOUBLE PRECISION NOT NULL,
    somaticNonsenseOrFrameshift DOUBLE PRECISION NOT NULL,
    somaticMissense DOUBLE PRECISION NOT NULL,
    somaticSplice DOUBLE PRECISION NOT NULL,
    somaticSynonymous DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (id),
    KEY(sampleId, type)
);

DROP TABLE IF EXISTS cuppaResult;
CREATE TABLE cuppaResult
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    cuppaResult varchar(255) NOT NULL,
    PRIMARY KEY (id),
    KEY(sampleId)
);


SET FOREIGN_KEY_CHECKS = 1;