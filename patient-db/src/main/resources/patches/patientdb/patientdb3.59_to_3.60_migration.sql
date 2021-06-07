DROP TABLE pgxCalls;

ALTER TABLE pgxGenotype RENAME TO peachGenotype;

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

CREATE TABLE cuppa
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    cuppaTumorLocation varchar(255) NOT NULL,
    cuppaPrediction varchar(255) NOT NULL,
    PRIMARY KEY (id),
    KEY(sampleId)
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

CREATE TABLE snpcheck
(   sampleId varchar(255) NOT NULL,
    isPass BOOLEAN NOT NULL,
    PRIMARY KEY (sampleId)
);

ALTER TABLE sample
    ADD COLUMN reportedDate DATE AFTER samplingDate;

DROP TABLE germlineVariant2;
DROP TABLE germlineVariant;

CREATE TABLE germlineVariant
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