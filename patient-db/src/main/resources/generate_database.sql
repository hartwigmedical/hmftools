-- noinspection SqlNoDataSourceInspectionForFile
-- All database and field names are wrapped in backticks to force jOOQ code generation to use them literally

DROP TABLE IF EXISTS `snpcheck`;
CREATE TABLE `snpcheck`
(   `sampleId` VARCHAR(255) NOT NULL,
    `isPass` TINYINT(1) NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `amberPatient`;
CREATE TABLE `amberPatient`
(   `modified` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    `sampleId` VARCHAR(255) NOT NULL,
    `patientId` INT NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `amberAnonymous`;
CREATE TABLE `amberAnonymous`
(   `modified` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    `sampleId` VARCHAR(255) NOT NULL,
    `hmfSampleId` VARCHAR(255) NOT NULL,
    `deleted` TINYINT(1) NOT NULL,
    PRIMARY KEY (`sampleId`)
);

CREATE INDEX `amberAnonymous_hmfSampleId` ON `amberAnonymous` (`hmfSampleId`);


DROP TABLE IF EXISTS `amberMapping`;
CREATE TABLE `amberMapping`
(   `modified` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    `firstSampleId` VARCHAR(255) NOT NULL,
    `secondSampleId` VARCHAR(255) NOT NULL,
    `matches` INT NOT NULL,
    `sites` INT NOT NULL,
    `likelihood` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`firstSampleId`, `secondSampleId`)
);

DROP TABLE IF EXISTS `amberSample`;
CREATE TABLE `amberSample`
(   `modified` DATETIME NOT NULL DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
    `sampleId` VARCHAR(255) NOT NULL,
    `site1` TINYINT NOT NULL,
    `site2` TINYINT NOT NULL,
    `site3` TINYINT NOT NULL,
    `site4` TINYINT NOT NULL,
    `site5` TINYINT NOT NULL,
    `site6` TINYINT NOT NULL,
    `site7` TINYINT NOT NULL,
    `site8` TINYINT NOT NULL,
    `site9` TINYINT NOT NULL,
    `site10` TINYINT NOT NULL,
    `site11` TINYINT NOT NULL,
    `site12` TINYINT NOT NULL,
    `site13` TINYINT NOT NULL,
    `site14` TINYINT NOT NULL,
    `site15` TINYINT NOT NULL,
    `site16` TINYINT NOT NULL,
    `site17` TINYINT NOT NULL,
    `site18` TINYINT NOT NULL,
    `site19` TINYINT NOT NULL,
    `site20` TINYINT NOT NULL,
    `site21` TINYINT NOT NULL,
    `site22` TINYINT NOT NULL,
    `site23` TINYINT NOT NULL,
    `site24` TINYINT NOT NULL,
    `site25` TINYINT NOT NULL,
    `site26` TINYINT NOT NULL,
    `site27` TINYINT NOT NULL,
    `site28` TINYINT NOT NULL,
    `site29` TINYINT NOT NULL,
    `site30` TINYINT NOT NULL,
    `site31` TINYINT NOT NULL,
    `site32` TINYINT NOT NULL,
    `site33` TINYINT NOT NULL,
    `site34` TINYINT NOT NULL,
    `site35` TINYINT NOT NULL,
    `site36` TINYINT NOT NULL,
    `site37` TINYINT NOT NULL,
    `site38` TINYINT NOT NULL,
    `site39` TINYINT NOT NULL,
    `site40` TINYINT NOT NULL,
    `site41` TINYINT NOT NULL,
    `site42` TINYINT NOT NULL,
    `site43` TINYINT NOT NULL,
    `site44` TINYINT NOT NULL,
    `site45` TINYINT NOT NULL,
    `site46` TINYINT NOT NULL,
    `site47` TINYINT NOT NULL,
    `site48` TINYINT NOT NULL,
    `site49` TINYINT NOT NULL,
    `site50` TINYINT NOT NULL,
    `site51` TINYINT NOT NULL,
    `site52` TINYINT NOT NULL,
    `site53` TINYINT NOT NULL,
    `site54` TINYINT NOT NULL,
    `site55` TINYINT NOT NULL,
    `site56` TINYINT NOT NULL,
    `site57` TINYINT NOT NULL,
    `site58` TINYINT NOT NULL,
    `site59` TINYINT NOT NULL,
    `site60` TINYINT NOT NULL,
    `site61` TINYINT NOT NULL,
    `site62` TINYINT NOT NULL,
    `site63` TINYINT NOT NULL,
    `site64` TINYINT NOT NULL,
    `site65` TINYINT NOT NULL,
    `site66` TINYINT NOT NULL,
    `site67` TINYINT NOT NULL,
    `site68` TINYINT NOT NULL,
    `site69` TINYINT NOT NULL,
    `site70` TINYINT NOT NULL,
    `site71` TINYINT NOT NULL,
    `site72` TINYINT NOT NULL,
    `site73` TINYINT NOT NULL,
    `site74` TINYINT NOT NULL,
    `site75` TINYINT NOT NULL,
    `site76` TINYINT NOT NULL,
    `site77` TINYINT NOT NULL,
    `site78` TINYINT NOT NULL,
    `site79` TINYINT NOT NULL,
    `site80` TINYINT NOT NULL,
    `site81` TINYINT NOT NULL,
    `site82` TINYINT NOT NULL,
    `site83` TINYINT NOT NULL,
    `site84` TINYINT NOT NULL,
    `site85` TINYINT NOT NULL,
    `site86` TINYINT NOT NULL,
    `site87` TINYINT NOT NULL,
    `site88` TINYINT NOT NULL,
    `site89` TINYINT NOT NULL,
    `site90` TINYINT NOT NULL,
    `site91` TINYINT NOT NULL,
    `site92` TINYINT NOT NULL,
    `site93` TINYINT NOT NULL,
    `site94` TINYINT NOT NULL,
    `site95` TINYINT NOT NULL,
    `site96` TINYINT NOT NULL,
    `site97` TINYINT NOT NULL,
    `site98` TINYINT NOT NULL,
    `site99` TINYINT NOT NULL,
    `site100` TINYINT NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `canonicalTranscript`;
CREATE TABLE `canonicalTranscript`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `assembly` VARCHAR(255) NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `geneId` VARCHAR(255) NOT NULL,
    `chromosomeBand` VARCHAR(255) NOT NULL,
    `chromosome` VARCHAR(10) NOT NULL,
    `geneStart` INT NOT NULL,
    `geneEnd` INT NOT NULL,
    `transcriptId` VARCHAR(255) NOT NULL,
    `transcriptStart` INT NOT NULL,
    `transcriptEnd` INT NOT NULL,
    `exons` INT NOT NULL,
    `exonStart` INT NOT NULL,
    `exonEnd` INT NOT NULL,
    `exonBases` INT NOT NULL,
    `strand` VARCHAR(255) NOT NULL,
    `codingStart` INT NOT NULL,
    `codingEnd` INT NOT NULL,
    `codingBases` INT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `canonicalTranscript_gene` ON `canonicalTranscript` (`gene`);
CREATE INDEX `canonicalTranscript_transcriptId` ON `canonicalTranscript` (`transcriptId`);

DROP TABLE IF EXISTS `driverGenePanel`;
CREATE TABLE `driverGenePanel`
(   `modified` DATETIME NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `reportMissense` TINYINT(1) NOT NULL,
    `reportNonsense` TINYINT(1) NOT NULL,
    `reportSplice` TINYINT(1) NOT NULL,
    `reportDeletion` TINYINT(1) NOT NULL,
    `reportDisruption` TINYINT(1) NOT NULL,
    `reportAmplification` TINYINT(1) NOT NULL,
    `reportSomaticHotspot` TINYINT(1) NOT NULL,
    `reportGermlineVariant` VARCHAR(50) NOT NULL,
    `reportGermlineHotspot` VARCHAR(50) NOT NULL,
    `likelihoodType` VARCHAR(255) NOT NULL,
    `reportGermlineDisruption` VARCHAR(50) NOT NULL,
    `reportGermlineDeletion` VARCHAR(50) NOT NULL,
    `additionalReportedTranscripts` VARCHAR(255) NOT NULL,
    `reportPGX` TINYINT(1) NOT NULL,
    PRIMARY KEY (`gene`)
);

DROP TABLE IF EXISTS `metric`;
CREATE TABLE `metric`
(   `sampleId` VARCHAR(255) NOT NULL,
    `refMeanCoverage` DOUBLE PRECISION NOT NULL,
    `refSdCoverage` DOUBLE PRECISION NOT NULL,
    `refMedianCoverage` INT NOT NULL,
    `refMadCoverage` INT NOT NULL,
    `refPctExcAdapter` DOUBLE PRECISION,
    `refPctExcMapQ` DOUBLE PRECISION NOT NULL,
    `refPctExcDupe` DOUBLE PRECISION NOT NULL,
    `refPctExcUnpaired` DOUBLE PRECISION NOT NULL,
    `refPctExcBaseQ` DOUBLE PRECISION NOT NULL,
    `refPctExcOverlap` DOUBLE PRECISION NOT NULL,
    `refPctExcCapped` DOUBLE PRECISION NOT NULL,
    `refPctExcTotal` DOUBLE PRECISION NOT NULL,
    `refCoverage1xPercentage` DOUBLE PRECISION,
    `refCoverage10xPercentage` DOUBLE PRECISION NOT NULL,
    `refCoverage20xPercentage` DOUBLE PRECISION NOT NULL,
    `refCoverage30xPercentage` DOUBLE PRECISION NOT NULL,
    `refCoverage60xPercentage` DOUBLE PRECISION NOT NULL,
    `tumorMeanCoverage` DOUBLE PRECISION NOT NULL,
    `tumorSdCoverage` DOUBLE PRECISION NOT NULL,
    `tumorMedianCoverage` INT NOT NULL,
    `tumorMadCoverage` INT NOT NULL,
    `tumorPctExcAdapter` DOUBLE PRECISION,
    `tumorPctExcMapQ` DOUBLE PRECISION NOT NULL,
    `tumorPctExcDupe` DOUBLE PRECISION NOT NULL,
    `tumorPctExcUnpaired` DOUBLE PRECISION NOT NULL,
    `tumorPctExcBaseQ` DOUBLE PRECISION NOT NULL,
    `tumorPctExcOverlap` DOUBLE PRECISION NOT NULL,
    `tumorPctExcCapped` DOUBLE PRECISION NOT NULL,
    `tumorPctExcTotal` DOUBLE PRECISION NOT NULL,
    `tumorCoverage1xPercentage` DOUBLE PRECISION,
    `tumorCoverage10xPercentage` DOUBLE PRECISION NOT NULL,
    `tumorCoverage20xPercentage` DOUBLE PRECISION NOT NULL,
    `tumorCoverage30xPercentage` DOUBLE PRECISION NOT NULL,
    `tumorCoverage60xPercentage` DOUBLE PRECISION NOT NULL,
    `sufficientCoverage` TINYINT(1) NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `flagstat`;
CREATE TABLE `flagstat`
(   `sampleId` VARCHAR(255) NOT NULL,
    `refUniqueReadCount` BIGINT NOT NULL,
    `refSecondaryCount` BIGINT NOT NULL,
    `refSupplementaryCount` BIGINT NOT NULL,
    `refDuplicateProportion` DOUBLE PRECISION NOT NULL,
    `refMappedProportion` DOUBLE PRECISION NOT NULL,
    `refPairedInSequencingProportion` DOUBLE PRECISION NOT NULL,
    `refProperlyPairedProportion` DOUBLE PRECISION NOT NULL,
    `refWithItselfAndMateMappedProportion` DOUBLE PRECISION NOT NULL,
    `refSingletonProportion` DOUBLE PRECISION NOT NULL,
    `tumorUniqueReadCount` BIGINT NOT NULL,
    `tumorSecondaryCount` BIGINT NOT NULL,
    `tumorSupplementaryCount` BIGINT NOT NULL,
    `tumorDuplicateProportion` DOUBLE PRECISION NOT NULL,
    `tumorMappedProportion` DOUBLE PRECISION NOT NULL,
    `tumorPairedInSequencingProportion` DOUBLE PRECISION NOT NULL,
    `tumorProperlyPairedProportion` DOUBLE PRECISION NOT NULL,
    `tumorWithItselfAndMateMappedProportion` DOUBLE PRECISION NOT NULL,
    `tumorSingletonProportion` DOUBLE PRECISION NOT NULL,
    `passQC` TINYINT(1) NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `somaticVariant`;
CREATE TABLE `somaticVariant`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(31) NOT NULL,
    `position` INT NOT NULL,
    `filter` VARCHAR(255) NOT NULL,
    `type` VARCHAR(10) NOT NULL,
    `ref` VARCHAR(255) NOT NULL,
    `alt` VARCHAR(255) NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `genesAffected` INT NOT NULL,
    `canonicalEffect` VARCHAR(255) NOT NULL,
    `canonicalCodingEffect` VARCHAR(50) NOT NULL,
    `canonicalHgvsCodingImpact` VARCHAR(100) NOT NULL,
    `canonicalHgvsProteinImpact` VARCHAR(200) NOT NULL,
    `spliceRegion` TINYINT(1) NOT NULL,
    `otherTranscriptEffects` VARCHAR(255) NOT NULL,
    `worstCodingEffect` VARCHAR(50) NOT NULL,
    `microhomology` VARCHAR(255) NOT NULL,
    `repeatSequence` VARCHAR(255) NOT NULL,
    `repeatCount` INT NOT NULL,
    `alleleReadCount` INT NOT NULL,
    `totalReadCount` INT NOT NULL,
    `adjustedVaf` DOUBLE PRECISION NOT NULL,
    `variantCopyNumber` DOUBLE PRECISION NOT NULL,
    `copyNumber` DOUBLE PRECISION NOT NULL,
    `tier` VARCHAR(20) NOT NULL,
    `trinucleotideContext` VARCHAR(3) NOT NULL,
    `subclonalLikelihood` DOUBLE PRECISION NOT NULL,
    `biallelic` TINYINT(1) NOT NULL,
    `hotspot` VARCHAR(20) NOT NULL,
    `mappability` DOUBLE PRECISION NOT NULL,
    `germlineStatus` VARCHAR(50) NOT NULL,
    `minorAlleleCopyNumber` DOUBLE PRECISION NOT NULL,
    `recovered` TINYINT(1) NOT NULL,
    `kataegis` VARCHAR(20) NOT NULL,
    `referenceAlleleReadCount` INT,
    `referenceTotalReadCount` INT,
    `rnaAlleleReadCount` INT,
    `rnaTotalReadCount` INT,
    `localPhaseSet` VARCHAR(50),
    `qual` DOUBLE PRECISION NOT NULL,
    `reported` TINYINT(1) NOT NULL,
    `clinvarInfo` VARCHAR(255) NULL,
    `gnomadFrequency` DOUBLE PRECISION NULL,
    `somaticLikelihood` VARCHAR(10) NULL,
    PRIMARY KEY (`id`)
);

CREATE UNIQUE INDEX `somaticVariant_unique_constraint` ON `somaticVariant` (`sampleId`, `chromosome`, `position`, `ref`, `alt`);
CREATE INDEX `somaticVariant_sampleId` ON `somaticVariant` (`sampleId`);
CREATE INDEX `somaticVariant_filter` ON `somaticVariant` (`filter`);
CREATE INDEX `somaticVariant_type` ON `somaticVariant` (`type`);
CREATE INDEX `somaticVariant_gene` ON `somaticVariant` (`gene`);

DROP TABLE IF EXISTS `germlineVariant`;
CREATE TABLE `germlineVariant`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(255) NOT NULL,
    `position` INT NOT NULL,
    `filter` VARCHAR(255) NOT NULL,
    `type` VARCHAR(255) NOT NULL,
    `ref` VARCHAR(255) NOT NULL,
    `alt` VARCHAR(255) NOT NULL,

    -- ### SAGE
    `qual` DOUBLE PRECISION NOT NULL,
    `tier` VARCHAR(20) NOT NULL,
    `germlineGenotype` VARCHAR(255) NOT NULL,
    `germlineAlleleReadCount` INT,
    `germlineTotalReadCount` INT,
    `rnaAlleleReadCount` INT,
    `rnaTotalReadCount` INT,
    `tumorAlleleReadCount` INT NOT NULL,
    `tumorTotalReadCount` INT NOT NULL,
    `localPhaseSet` VARCHAR(50),

    -- ### PURPLE ENRICHMENT
    `adjustedVaf` DOUBLE PRECISION NOT NULL,
    `variantCopyNumber` DOUBLE PRECISION NOT NULL,
    `copyNumber` DOUBLE PRECISION NOT NULL,
    `biallelic` TINYINT(1) NOT NULL,
    `minorAlleleCopyNumber` DOUBLE PRECISION NOT NULL,

    -- ### PATHOGENIC
    `clinvarInfo` VARCHAR(255) NOT NULL,
    `pathogenicity` VARCHAR(255) NOT NULL,
    `pathogenic` TINYINT(1) NOT NULL,

    -- ### PAVE IMPACTS
    `gene` VARCHAR(50) NOT NULL,
    `genesAffected` INT NOT NULL,
    `canonicalEffect` VARCHAR(255) NOT NULL,
    `canonicalCodingEffect` VARCHAR(50) NOT NULL,
    `canonicalHgvsCodingImpact` VARCHAR(100) NOT NULL,
    `canonicalHgvsProteinImpact` VARCHAR(200) NOT NULL,
    `spliceRegion` TINYINT(1) NOT NULL,
    `otherTranscriptEffects` VARCHAR(255) NOT NULL,
    `worstCodingEffect` VARCHAR(50) NOT NULL,

    -- ### REF GENOME ENRICHMENT
    `microhomology` VARCHAR(255) NOT NULL,
    `repeatSequence` VARCHAR(255) NOT NULL,
    `repeatCount` INT NOT NULL,
    `trinucleotideContext` VARCHAR(3) NOT NULL,

    -- ### DRIVER CATALOG
    `hotspot` VARCHAR(20) NOT NULL,
    `mappability` DOUBLE PRECISION NOT NULL,
    `reported` TINYINT(1) NOT NULL,

    PRIMARY KEY (`id`)
);

CREATE INDEX `germlineVariant_sampleId` ON `germlineVariant` (`sampleId`);
CREATE INDEX `germlineVariant_filter` ON `germlineVariant` (`filter`);
CREATE INDEX `germlineVariant_type` ON `germlineVariant` (`type`);
CREATE INDEX `germlineVariant_gene` ON `germlineVariant` (`gene`);

DROP TABLE IF EXISTS `purity`;
CREATE TABLE `purity`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `gender` VARCHAR(20) NOT NULL,
    `runMode` VARCHAR(20) NOT NULL,
    `fitMethod` VARCHAR(20) NOT NULL,
    `qcStatus` VARCHAR(255) NOT NULL,
    `purity` DOUBLE PRECISION NOT NULL,
    `normFactor` DOUBLE PRECISION NOT NULL,
    `score` DOUBLE PRECISION NOT NULL,
    `somaticPenalty` DOUBLE PRECISION NOT NULL,
    `ploidy` DOUBLE PRECISION NOT NULL,
    `diploidProportion` DOUBLE PRECISION NOT NULL,
    `polyclonalProportion` DOUBLE PRECISION NOT NULL,
    `wholeGenomeDuplication` TINYINT(1) NOT NULL,
    `minPurity` DOUBLE PRECISION NOT NULL,
    `maxPurity` DOUBLE PRECISION NOT NULL,
    `minPloidy` DOUBLE PRECISION NOT NULL,
    `maxPloidy` DOUBLE PRECISION NOT NULL,
    `minDiploidProportion` DOUBLE PRECISION NOT NULL,
    `maxDiploidProportion` DOUBLE PRECISION NOT NULL,
    `msIndelsPerMb` DOUBLE PRECISION NOT NULL,
    `msStatus` VARCHAR(10) NOT NULL,
    `tmbPerMb` DOUBLE PRECISION NOT NULL,
    `tmbStatus` VARCHAR(10) NOT NULL,
    `tml` INT NOT NULL,
    `tmlStatus` VARCHAR(10) NOT NULL,
    `svTmb` INT NOT NULL DEFAULT 0,
    `deletedGenes` INT NOT NULL DEFAULT 0,
    `copyNumberSegments` INT NOT NULL DEFAULT 0,
    `unsupportedCopyNumberSegments` INT NOT NULL DEFAULT 0,
    `contamination` DOUBLE PRECISION NOT NULL DEFAULT 0,
    `germlineAberration` VARCHAR(255) NOT NULL DEFAULT 'NONE',
    `amberGender` VARCHAR(20) NOT NULL,
    `targeted` TINYINT(1) NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `purity_sampleId` ON `purity` (`sampleId`);

DROP TABLE IF EXISTS `purityRange`;
CREATE TABLE `purityRange`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `purity` DOUBLE PRECISION NOT NULL,
    `normFactor` DOUBLE PRECISION NOT NULL,
    `score` DOUBLE PRECISION NOT NULL,
    `somaticPenalty` DOUBLE PRECISION NOT NULL,
    `ploidy` DOUBLE PRECISION NOT NULL,
    `diploidProportion` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `purityRange_sampleId` ON `purityRange` (`sampleId`);

DROP TABLE IF EXISTS `copyNumber`;
CREATE TABLE `copyNumber`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(255) NOT NULL,
    `start` INT NOT NULL,
    `end` INT NOT NULL,
    `segmentStartSupport` VARCHAR(255) NOT NULL,
    `segmentEndSupport` VARCHAR(255) NOT NULL,
    `depthWindowCount` INT NOT NULL,
    `bafCount` INT NOT NULL,
    `observedBaf` DOUBLE PRECISION NOT NULL,
    `baf` DOUBLE PRECISION NOT NULL,
    `copyNumber` DOUBLE PRECISION NOT NULL,
    `minorAlleleCopyNumber` DOUBLE PRECISION NOT NULL,
    `majorAlleleCopyNumber` DOUBLE PRECISION NOT NULL,
    `copyNumberMethod` VARCHAR(255) NOT NULL,
    `gcContent` DOUBLE PRECISION NOT NULL,
    `minStart` INT NOT NULL,
    `maxStart` INT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE UNIQUE INDEX `copyNumber_unique_constraint` ON `copyNumber` (`sampleId`, `chromosome`, `start`, `end`);
CREATE INDEX `copyNumber_sampleId` ON `copyNumber` (`sampleId`);

DROP TABLE IF EXISTS `geneCopyNumber`;
CREATE TABLE `geneCopyNumber`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(31) NOT NULL,
    `start` INT NOT NULL,
    `end` INT NOT NULL,
    `gene` VARCHAR(255) NOT NULL,
    `chromosomeBand` VARCHAR(255) NOT NULL,
    `transcriptId` VARCHAR(255) NOT NULL,
    `canonicalTranscript` TINYINT(1) NOT NULL,
    `minCopyNumber` DOUBLE PRECISION NOT NULL,
    `maxCopyNumber` DOUBLE PRECISION NOT NULL,
    `somaticRegions` INT NOT NULL,
    `minRegions` INT NOT NULL,
    `minRegionStart` INT NOT NULL,
    `minRegionEnd` INT NOT NULL,
    `minRegionStartSupport` VARCHAR(255) NOT NULL,
    `minRegionEndSupport` VARCHAR(255) NOT NULL,
    `minRegionMethod` VARCHAR(255) NOT NULL,
    `minMinorAlleleCopyNumber` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE UNIQUE INDEX `geneCopyNumber_unique_constraint` ON `geneCopyNumber` (`sampleId`, `chromosome`, `gene`, `transcriptId`);
CREATE INDEX `geneCopyNumber_sampleId_gene` ON `geneCopyNumber` (`sampleId`, `gene`);
CREATE INDEX `geneCopyNumber_gene` ON `geneCopyNumber` (`gene`);

DROP TABLE IF EXISTS `driverCatalog`;
CREATE TABLE `driverCatalog`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(10) NOT NULL,
    `chromosomeBand` VARCHAR(50) NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `transcriptId` VARCHAR(50) NOT NULL,
    `canonicalTranscript` TINYINT(1) NOT NULL DEFAULT 1,
    `category` VARCHAR(50) NOT NULL,
    `driver` VARCHAR(50) NOT NULL,
    `likelihoodMethod` VARCHAR(50) NOT NULL,
    `driverLikelihood` DOUBLE PRECISION NOT NULL,
    `missense` INT NOT NULL,
    `nonsense` INT NOT NULL,
    `splice` INT NOT NULL,
    `frameshift` INT NOT NULL,
    `inframe` INT NOT NULL,
    `biallelic` TINYINT(1) NOT NULL,
    `minCopyNumber` DOUBLE PRECISION NOT NULL,
    `maxCopyNumber` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `driverCatalog_sampleId_gene` ON `driverCatalog` (`sampleId`, `gene`);
CREATE INDEX `driverCatalog_gene` ON `driverCatalog` (`gene`);

DROP TABLE IF EXISTS `germlineDeletion`;
CREATE TABLE `germlineDeletion`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `chromosome` VARCHAR(10) NOT NULL,
    `chromosomeBand` VARCHAR(50) NOT NULL,
    `regionStart` INT NOT NULL,
    `regionEnd` INT NOT NULL,
    `depthWindowCount` INT NOT NULL,
    `exonStart` INT NOT NULL,
    `exonEnd` INT NOT NULL,
    `detectionMethod` VARCHAR(20) NOT NULL,
    `germlineStatus` VARCHAR(20) NOT NULL,
    `tumorStatus` VARCHAR(20) NOT NULL,
    `germlineCopyNumber` DOUBLE PRECISION NOT NULL,
    `tumorCopyNumber` DOUBLE PRECISION NOT NULL,
    `filter` VARCHAR(50) NOT NULL,
    `cohortFrequency` INT NOT NULL,
    `reported` TINYINT(1) NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `germlineDeletion_sampleId_gene` ON `germlineDeletion` (`sampleId`, `gene`);
CREATE INDEX `germlineDeletion_gene` ON `germlineDeletion` (`gene`);

DROP TABLE IF EXISTS `structuralVariant`;
CREATE TABLE `structuralVariant`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `svId` INT NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `startChromosome` VARCHAR(31) NOT NULL,
    `endChromosome` VARCHAR(31),
    `startPosition` INT NOT NULL,
    `endPosition` INT,
    `startOrientation` TINYINT NOT NULL,
    `endOrientation` TINYINT,
    `startHomologySequence` VARCHAR(255) NOT NULL,
    `endHomologySequence` VARCHAR(255),
    `startAF` DOUBLE PRECISION,
    `endAF` DOUBLE PRECISION,
    `junctionCopyNumber` DOUBLE PRECISION,
    `adjustedAFStart` DOUBLE PRECISION,
    `adjustedAFEnd` DOUBLE PRECISION,
    `adjustedCopyNumberStart` DOUBLE PRECISION,
    `adjustedCopyNumberEnd` DOUBLE PRECISION,
    `adjustedCopyNumberChangeStart` DOUBLE PRECISION,
    `adjustedCopyNumberChangeEnd` DOUBLE PRECISION,
    `insertSequence` VARCHAR(2048) NOT NULL,
    `type` VARCHAR(10) NOT NULL,
    `filter` VARCHAR(50) NOT NULL,
    `imprecise` TINYINT(1) NOT NULL,
    `qualScore` DOUBLE PRECISION,
    `event` VARCHAR(50),
    `startTumorVariantFragmentCount` INT,
    `startTumorReferenceFragmentCount` INT,
    `startNormalVariantFragmentCount` INT,
    `startNormalReferenceFragmentCount` INT,
    `endTumorVariantFragmentCount` INT,
    `endTumorReferenceFragmentCount` INT,
    `endNormalVariantFragmentCount` INT,
    `endNormalReferenceFragmentCount` INT,
    `startIntervalOffsetStart` INT,
    `startIntervalOffsetEnd` INT,
    `endIntervalOffsetStart` INT,
    `endIntervalOffsetEnd` INT,
    `inexactHomologyOffsetStart` INT,
    `inexactHomologyOffsetEnd` INT,
    `startLinkedBy` VARCHAR(1024),
    `endLinkedBy` VARCHAR(1024),
    `vcfId` VARCHAR(50),
    `recovered` TINYINT(1) NOT NULL,
    `recoveryMethod` VARCHAR(64),
    `recoveryFilter` VARCHAR(255),
    `startRefContext` VARCHAR(255),
    `endRefContext` VARCHAR(255),
    `insertSequenceAlignments` VARCHAR(512),
    `insertSequenceRepeatClass` VARCHAR(64),
    `insertSequenceRepeatType` VARCHAR(64),
    `insertSequenceRepeatOrientation` TINYINT,
    `insertSequenceRepeatCoverage` DOUBLE PRECISION,
    `startAnchoringSupportDistance` INT,
    `endAnchoringSupportDistance` INT,
    PRIMARY KEY (`id`)
);

CREATE INDEX `structuralVariant_sampleId_svId` ON `structuralVariant` (`sampleId`, `svId`);

DROP TABLE IF EXISTS `svAnnotation`;
CREATE TABLE `svAnnotation`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `svId` INT NOT NULL,
    `clusterId` INT NOT NULL,
    `clusterReason` VARCHAR(255) NULL,
    `fragileSiteStart` TINYINT(1) NOT NULL,
    `fragileSiteEnd` TINYINT(1) NOT NULL,
    `isFoldback` TINYINT(1) NOT NULL,
    `lineTypeStart` VARCHAR(20),
    `lineTypeEnd` VARCHAR(20),
    `junctionCopyNumberMin` DOUBLE PRECISION,
    `junctionCopyNumberMax` DOUBLE PRECISION,
    `geneStart` VARCHAR(20),
    `geneEnd` VARCHAR(20),
    `localTopologyIdStart` INT,
    `localTopologyIdEnd` INT,
    `localTopologyStart` VARCHAR(20),
    `localTopologyEnd` VARCHAR(20),
    `localTICountStart` INT,
    `localTICountEnd` INT,
    PRIMARY KEY (`id`)
);

CREATE INDEX `svAnnotation_sampleId_svId` ON `svAnnotation` (`sampleId`, `svId`);
CREATE INDEX `svAnnotation_clusterId` ON `svAnnotation` (`clusterId`);
CREATE INDEX `svAnnotation_svId` ON `svAnnotation` (`svId`);

DROP TABLE IF EXISTS `svCluster`;
CREATE TABLE `svCluster`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `clusterId` INT NOT NULL,
    `category` VARCHAR(20),
    `synthetic` TINYINT(1) NOT NULL,
    `resolvedType` VARCHAR(20),
    `clusterCount` INT,
    `clusterDesc` VARCHAR(50),
    PRIMARY KEY (`id`)
);

CREATE INDEX `svCluster_sampleId_clusterId` ON `svCluster` (`sampleId`, `clusterId`);
CREATE INDEX `svCluster_clusterId` ON `svCluster` (`clusterId`);

DROP TABLE IF EXISTS `svLink`;
CREATE TABLE `svLink`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `clusterId` INT NOT NULL,
    `chainId` INT NOT NULL,
    `chainIndex` VARCHAR(50) NOT NULL,
    `chainLinkCount` INT NOT NULL,
    `lowerSvId` INT NOT NULL,
    `upperSvId` INT NOT NULL,
    `lowerBreakendIsStart` TINYINT(1) NOT NULL,
    `upperBreakendIsStart` TINYINT(1) NOT NULL,
    `chromosome` VARCHAR(10),
    `arm` VARCHAR(3),
    `assembled` TINYINT(1) NOT NULL,
    `traversedSVCount` INT,
    `linkLength` INT,
    `junctionCopyNumber` DOUBLE PRECISION,
    `junctionCopyNumberUncertainty` DOUBLE PRECISION,
    `pseudogeneInfo` VARCHAR(255),
    `ecDna` TINYINT(1),
    PRIMARY KEY (`id`)
);

CREATE INDEX `svLink_sampleId_clusterId` ON `svLink` (`sampleId`, `clusterId`);
CREATE INDEX `svLink_clusterId` ON `svLink` (`clusterId`);

DROP TABLE IF EXISTS `svDriver`;
CREATE TABLE `svDriver`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `clusterId` INT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `eventType` VARCHAR(50),
    PRIMARY KEY (`id`)
);

CREATE INDEX `svDiver_sampleId_clusterId` ON `svDriver` (`sampleId`, `clusterId`);
CREATE INDEX `svDriver_clusterId` ON `svDriver` (`clusterId`);

DROP TABLE IF EXISTS `svBreakend`;
CREATE TABLE `svBreakend`
(   `id` INT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `svId` INT NOT NULL,
    `startBreakend` TINYINT(1) NOT NULL,
    `gene` VARCHAR(50) NOT NULL, -- # length here comes from ensembl db schema
    `transcriptId` VARCHAR(128) NOT NULL, -- # length here comes from ensembl db schema
    `canonicalTranscript` TINYINT(1) NOT NULL,
    `geneOrientation` VARCHAR(20) NOT NULL,
    `disruptive` TINYINT(1) NOT NULL,
    `reportedDisruption` TINYINT(1) NOT NULL,
    `undisruptedCopyNumber` DOUBLE PRECISION,
    `regionType` VARCHAR(20) NOT NULL,
    `codingContext` VARCHAR(20),
    `biotype` VARCHAR(255),
    `exonUp` SMALLINT NOT NULL,
    `exonDown` SMALLINT NOT NULL,
    `exonicBasePhase` TINYINT,
    `nextSpliceExonRank` SMALLINT,
    `nextSpliceExonPhase` TINYINT,
    `nextSpliceDistance` INT,
    `totalExonCount` SMALLINT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `svBreakend_sampleId_svId` ON `svBreakend` (`sampleId`, `svId`);
CREATE INDEX `svBreakend_gene` ON `svBreakend` (`gene`);
CREATE INDEX `svBreakend_transcriptId` ON `svBreakend` (`transcriptId`);

DROP TABLE IF EXISTS `svFusion`;
CREATE TABLE `svFusion`
(   `id` INT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `fivePrimeBreakendId` INT UNSIGNED NOT NULL,
    `threePrimeBreakendId` INT UNSIGNED NOT NULL,
    `name` VARCHAR(50) NOT NULL,
    `reported` TINYINT(1) NOT NULL,
    `reportedType` VARCHAR(50) NOT NULL,
    `reportedReason` VARCHAR(255),
    `phased` VARCHAR(20) NOT NULL,
    `likelihood` VARCHAR(10) NOT NULL,
    `chainLength` INT,
    `chainLinks` INT,
    `chainTerminated` TINYINT(1),
    `domainsKept` VARCHAR(255),
    `domainsLost` VARCHAR(255),
    `skippedExonsUp` INT,
    `skippedExonsDown` INT,
    `fusedExonUp` INT,
    `fusedExonDown` INT,
    PRIMARY KEY (`id`)
);

CREATE INDEX `svFusion_fivePrimeBreakendId` ON `svFusion` (`fivePrimeBreakendId`);
CREATE INDEX `svFusion_threePrimeBreakendId` ON `svFusion` (`threePrimeBreakendId`);
CREATE INDEX `svFusion_sampleId` ON `svFusion` (`sampleId`);

DROP TABLE IF EXISTS `structuralVariantGermline`;
CREATE TABLE `structuralVariantGermline`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `svId` INT NOT NULL,
    `vcfId` VARCHAR(50),
    `chromosomeStart` VARCHAR(10) NOT NULL,
    `chromosomeEnd` VARCHAR(10),
    `positionStart` INT NOT NULL,
    `positionEnd` INT,
    `orientationStart` TINYINT NOT NULL,
    `orientationEnd` TINYINT,
    `type` VARCHAR(10) NOT NULL,
    `filter` VARCHAR(50) NOT NULL,
    `event` VARCHAR(50),
    `qualScore` DOUBLE PRECISION,
    `homologySequenceStart` VARCHAR(255) NOT NULL,
    `homologySequenceEnd` VARCHAR(255),
    `junctionCopyNumber` DOUBLE PRECISION,
    `adjustedAFStart` DOUBLE PRECISION,
    `adjustedAFEnd` DOUBLE PRECISION,
    `adjustedCopyNumberStart` DOUBLE PRECISION,
    `adjustedCopyNumberEnd` DOUBLE PRECISION,
    `adjustedCopyNumberChangeStart` DOUBLE PRECISION,
    `adjustedCopyNumberChangeEnd` DOUBLE PRECISION,
    `germlineFragments` INT,
    `germlineReferenceFragmentsStart` INT,
    `germlineReferenceFragmentsEnd` INT,
    `tumorFragments` INT,
    `tumorReferenceFragmentsStart` INT,
    `tumorReferenceFragmentsEnd` INT,
    `insertSequence` VARCHAR(2048) NOT NULL,
    `insertSequenceAlignments` VARCHAR(512),
    `insertSequenceRepeatClass` VARCHAR(64),
    `insertSequenceRepeatType` VARCHAR(64),
    `clusterId` INT NULL,
    `clusterCount` INT NULL,
    `resolvedType` VARCHAR(20),
    `linkedByStart` VARCHAR(1024),
    `linkedByEnd` VARCHAR(1024),
    `cohortFrequency` INT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `structuralVariantGermline_sampleId` ON `structuralVariantGermline` (`sampleId`);

DROP TABLE IF EXISTS `svBreakendGermline`;
CREATE TABLE `svBreakendGermline`
(   `id` INT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `svId` INT NOT NULL,
    `startBreakend` TINYINT(1) NOT NULL,
    `gene` VARCHAR(50) NOT NULL, -- length here comes from ensembl db schema
    `transcriptId` VARCHAR(128) NOT NULL, -- length here comes from ensembl db schema
    `canonicalTranscript` TINYINT(1) NOT NULL,
    `geneOrientation` VARCHAR(20) NOT NULL,
    `disruptive` TINYINT(1) NOT NULL,
    `reportedDisruption` TINYINT(1) NOT NULL,
    `undisruptedCopyNumber` DOUBLE PRECISION,
    `regionType` VARCHAR(20) NOT NULL,
    `codingType` VARCHAR(20),
    `biotype` VARCHAR(255),
    `exonicBasePhase` TINYINT,
    `nextSpliceExonRank` TINYINT UNSIGNED,
    `nextSpliceExonPhase` TINYINT,
    `nextSpliceDistance` INT,
    `totalExonCount` SMALLINT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `svBreakendGermline_sampleId_svId` ON `svBreakendGermline` (`sampleId`, `svId`);
CREATE INDEX `svBreakendGermline_gene` ON `svBreakendGermline` (`gene`);
CREATE INDEX `svBreakendGermline_transcriptId` ON `svBreakendGermline` (`transcriptId`);

DROP TABLE IF EXISTS `signature`;
CREATE TABLE `signature`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `signature` VARCHAR(50) NOT NULL,
    `allocation` DOUBLE PRECISION,
    `percent` DOUBLE PRECISION,
    PRIMARY KEY (`id`)
);

CREATE INDEX `signature_sampleId` ON `signature` (`sampleId`);

DROP TABLE IF EXISTS `chord`;
CREATE TABLE `chord`
(   `sampleId` VARCHAR(255) NOT NULL,
    `BRCA1` DOUBLE PRECISION NOT NULL,
    `BRCA2` DOUBLE PRECISION NOT NULL,
    `hrd` DOUBLE PRECISION NOT NULL,
    `hrStatus` VARCHAR(255) NOT NULL,
    `hrdType` VARCHAR(255) NOT NULL,
    `remarksHrStatus` VARCHAR(255) NOT NULL,
    `remarksHrdType` VARCHAR(255) NOT NULL,
    PRIMARY KEY (`sampleId`)
);

DROP TABLE IF EXISTS `peachGenotype`;
CREATE TABLE `peachGenotype`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(255) NOT NULL,
    `gene` VARCHAR(255) NOT NULL,
    `haplotype` VARCHAR(255) NOT NULL,
    `count` INT NOT NULL,
    `function` VARCHAR(255) NOT NULL,
    `linkedDrugs` VARCHAR(255) NOT NULL,
    `urlPrescriptionInfo` VARCHAR(255) NOT NULL,
    PRIMARY KEY (`id`)
);

DROP TABLE IF EXISTS `virusBreakend`;
CREATE TABLE `virusBreakend`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(255) NOT NULL,
    `taxidGenus` INT NOT NULL,
    `nameGenus` VARCHAR(255) NOT NULL,
    `readsGenusTree` INT NOT NULL,
    `taxidSpecies` INT NOT NULL,
    `nameSpecies` VARCHAR(255) NOT NULL,
    `readsSpeciesTree` INT NOT NULL,
    `taxidAssigned` INT NOT NULL,
    `nameAssigned` VARCHAR(255) NOT NULL,
    `readsAssignedTree` INT NOT NULL,
    `readsAssignedDirect` INT NOT NULL,
    `reference` VARCHAR(255) NOT NULL,
    `referenceTaxid` INT NOT NULL,
    `referenceKmerCount` INT NOT NULL,
    `alternateKmerCount` INT NOT NULL,
    `RName` VARCHAR(255) NOT NULL,
    `startPos` INT NOT NULL,
    `endPos` INT NOT NULL,
    `numReads` INT NOT NULL,
    `covBases` INT NOT NULL,
    `coverage` INT NOT NULL,
    `meanDepth` INT NOT NULL,
    `meanBaseQ` INT NOT NULL,
    `meanMapQ` INT NOT NULL,
    `integrations` INT NOT NULL,
    `qcStatus` VARCHAR(255) NOT NULL,
    PRIMARY KEY (`id`)
);

DROP TABLE IF EXISTS `virusAnnotation`;
CREATE TABLE `virusAnnotation`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(255) NOT NULL,
    `taxid` INT NOT NULL,
    `virusName` VARCHAR(255) NOT NULL,
    `qcStatus` VARCHAR(255) NOT NULL,
    `integrations` INT NOT NULL,
    `interpretation` VARCHAR(255),
    `percentageCovered` DOUBLE NOT NULL,
    `meanCoverage` DOUBLE NOT NULL,
    `expectedClonalCoverage` DOUBLE,
    `reported` TINYINT(1) NOT NULL,
    `likelihood` VARCHAR(10) NOT NULL,
    PRIMARY KEY (`id`)
);

DROP TABLE IF EXISTS `cuppa`;
CREATE TABLE `cuppa`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(255) NOT NULL,
    `clfName` VARCHAR(255) DEFAULT NULL,
    `cancerType` VARCHAR(255),
    `prob` DOUBLE PRECISION,
    `rank` INT NOT NULL,
    `isOldCuppaOutput` TINYINT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `cuppa_sampleId` ON `cuppa` (`sampleId`);

DROP TABLE IF EXISTS `rnaStatistics`;
CREATE TABLE `rnaStatistics`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `qcStatus` VARCHAR(100) NOT NULL,
    `readLength` INT NOT NULL,
    `totalFragments` INT NOT NULL,
    `duplicates` INT NOT NULL,
    `splicedPercent` DOUBLE PRECISION NOT NULL,
    `unsplicedPercent` DOUBLE PRECISION NOT NULL,
    `alternateSplicePercent` DOUBLE PRECISION NOT NULL,
    `chimericPercent` DOUBLE PRECISION NOT NULL,
    `fragmentLengthPct05` DOUBLE PRECISION NOT NULL,
    `fragmentLengthPct50` DOUBLE PRECISION NOT NULL,
    `fragmentLengthPct95` DOUBLE PRECISION NOT NULL,
    `enrichedGenePercent` DOUBLE PRECISION NOT NULL,
    `medianGCRatio` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `rnaStatistics_sampleId` ON `rnaStatistics` (`sampleId`);

DROP TABLE IF EXISTS `geneExpression`;
CREATE TABLE `geneExpression`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `gene` VARCHAR(30) NOT NULL,
    `tpm` DOUBLE PRECISION NOT NULL,
    `splicedFragments` INT NOT NULL,
    `unsplicedFragments` INT NOT NULL,
    `medianTpmCancer` DOUBLE PRECISION NOT NULL,
    `percentileCancer` DOUBLE PRECISION NOT NULL,
    `medianTpmCohort` DOUBLE PRECISION NOT NULL,
    `percentileCohort` DOUBLE PRECISION NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `geneExpression_sampleId_gene` ON `geneExpression` (`sampleId`, `gene`);
CREATE INDEX `geneExpression_gene` ON `geneExpression` (`gene`);

DROP TABLE IF EXISTS `novelSpliceJunction`;
CREATE TABLE `novelSpliceJunction`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `gene` VARCHAR(20) NOT NULL,
    `chromosome` VARCHAR(10) NOT NULL,
    `junctionStart` INT NOT NULL,
    `junctionEnd` INT NOT NULL,
    `type` VARCHAR(20) NOT NULL,
    `fragmentCount` INT NOT NULL,
    `depthStart` INT NOT NULL,
    `depthEnd` INT NOT NULL,
    `regionStart` VARCHAR(20) NOT NULL,
    `regionEnd` VARCHAR(20) NOT NULL,
    `basesStart` VARCHAR(20) NOT NULL,
    `basesEnd` VARCHAR(20) NOT NULL,
    `cohortFrequency` INT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `novelSpliceJunction_sampleId_gene` ON `novelSpliceJunction` (`sampleId`, `gene`);
CREATE INDEX `novelSpliceJunction_gene` ON `novelSpliceJunction` (`gene`);

DROP TABLE IF EXISTS `rnaFusion`;
CREATE TABLE `rnaFusion`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `name` VARCHAR(50) NOT NULL,
    `chromosomeUp` VARCHAR(10) NOT NULL,
    `chromosomeDown` VARCHAR(10) NOT NULL,
    `positionUp` INT NOT NULL,
    `positionDown` INT NOT NULL,
    `orientationUp` TINYINT NOT NULL,
    `orientationDown` TINYINT NOT NULL,
    `junctionTypeUp` VARCHAR(20) NOT NULL,
    `junctionTypeDown` VARCHAR(20) NOT NULL,
    `svType` VARCHAR(5) NOT NULL,
    `splitFragments` INT NOT NULL,
    `realignedFrags` INT NOT NULL,
    `discordantFrags` INT NOT NULL,
    `depthUp` INT NOT NULL,
    `depthDown` INT NOT NULL,
    `maxAnchorLengthUp` INT NOT NULL,
    `maxAnchorLengthDown` INT NOT NULL,
    `cohortFrequency` INT NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `rnaFusion_sampleId_name` ON `rnaFusion` (`sampleId`, `name`);
CREATE INDEX `rnaFusion_name` ON `rnaFusion` (`name`);

DROP TABLE IF EXISTS `cdr3Sequence`;
CREATE TABLE `cdr3Sequence`
(   `id` INT NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `cdr3Seq` VARCHAR(255) NOT NULL,
    `cdr3AA` VARCHAR(100) NOT NULL,
    `locus` VARCHAR(10) NOT NULL,
    `filter` VARCHAR(50) NOT NULL,
    `blastnStatus` VARCHAR(10) NOT NULL,
    `minHighQualBaseReads` INT NOT NULL,
    `assignedReads` INT NOT NULL,
    `inFrame` TINYINT NOT NULL,
    `containsStop` TINYINT NOT NULL,
    PRIMARY KEY (`id`)
);

DROP TABLE IF EXISTS `cdr3LocusSummary`;
CREATE TABLE `cdr3LocusSummary`
(    `id` INT NOT NULL AUTO_INCREMENT,
     `modified` DATETIME NOT NULL,
     `sampleId` VARCHAR(50) NOT NULL,
     `locus` VARCHAR(10) NOT NULL,
     `readsUsed` INT NOT NULL,
     `readsTotal` INT NOT NULL,
     `downSampled` TINYINT NOT NULL,
     `sequences` INT NOT NULL,
     `passSequences` INT NOT NULL,
     PRIMARY KEY (`id`)
);

CREATE UNIQUE INDEX `cdr3LocusSummary_sampleId_locus` ON `cdr3LocusSummary` (`sampleId`, `locus`);

DROP TABLE IF EXISTS `telomereLength`;
CREATE TABLE `telomereLength`
(    `id` INT NOT NULL AUTO_INCREMENT,
     `modified` DATETIME NOT NULL,
     `sampleId` VARCHAR(50) NOT NULL,
     `germlineTelomereLength` DOUBLE PRECISION,
     `somaticTelomereLength` DOUBLE PRECISION,
     `germlineFullFragments` INT,
     `germlineCRichPartialFragments` INT,
     `germlineGRichPartialFragments` INT,
     `somaticFullFragments` INT,
     `somaticCRichPartialFragments` INT,
     `somaticGRichPartialFragments` INT,
     `sampleMixLength` DOUBLE PRECISION,
     PRIMARY KEY (`id`)
);

CREATE UNIQUE INDEX `telomereLength_sampleId` ON `telomereLength` (`sampleId`);
