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
    Rname varchar(255) NOT NULL,
    startpos int NOT NULL,
    endpos int NOT NULL,
    numreads int NOT NULL,
    covbases int NOT NULL,
    coverage int NOT NULL,
    meandepth int NOT NULL,
    meanbaseq int NOT NULL,
    meanmapq int NOT NULL,
    integrations int NOT NULL,
    QCStatus varchar(255) NOT NULL,
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

CREATE TABLE cuppaResult
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    sampleId varchar(255) NOT NULL,
    cuppaResult varchar(255) NOT NULL,
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