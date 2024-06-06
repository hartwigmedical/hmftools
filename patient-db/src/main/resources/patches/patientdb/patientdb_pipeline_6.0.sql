####
# SQL updates for Pipeline release 6.0
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# DEV-2919: Java rewrite Peach

DROP TABLE IF EXISTS peachCalls;

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



# Linx

ALTER TABLE somaticVariant
    ADD COLUMN reportedReason varchar(255) NULL after reportedType;
