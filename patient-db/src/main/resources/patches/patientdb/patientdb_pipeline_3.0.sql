####
# SQL updates for Pipeline release 3.0

# Drivers
ALTER TABLE driverCatalog
    ADD COLUMN reportedStatus varchar(50) NULL after likelihoodMethod;

ALTER TABLE svBreakend
    DROP COLUMN reportedDisruption,
    ADD COLUMN reportedStatus varchar(50) NULL after disruptive;

ALTER TABLE svBreakend
    DROP COLUMN svBreakendGermline,
    ADD COLUMN reportedStatus varchar(50) NULL after disruptive;

ALTER TABLE somaticVariant
    DROP COLUMN recovered;

# or rename to keep records
DROP TABLE IF EXISTS `germlineCopyNumber`;
CREATE TABLE `germlineCopyNumber`
(   `id` BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
    `modified` DATETIME NOT NULL,
    `sampleId` VARCHAR(50) NOT NULL,
    `gene` VARCHAR(50) NOT NULL,
    `transcriptId` VARCHAR(128) NOT NULL,
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
    `normalCopyNumber` DOUBLE PRECISION NOT NULL,
    `tumorCopyNumber` DOUBLE PRECISION NOT NULL,
    `filter` VARCHAR(50) NOT NULL,
    `cohortFrequency` INT NOT NULL,
    `reported` TINYINT(1) NOT NULL,
    PRIMARY KEY (`id`)
);

CREATE INDEX `germlineCopyNumber_sampleId_gene` ON `germlineCopyNumber` (`sampleId`, `gene`);
CREATE INDEX `germlineCopyNumber_gene` ON `germlineCopyNumber` (`gene`);