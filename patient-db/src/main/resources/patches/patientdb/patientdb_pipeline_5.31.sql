####
# SQL updates for Pipeline release 5.31
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Protect

ALTER TABLE protect
    ADD COLUMN sourceTreatmentApproach varchar(500) after treatment;

ALTER TABLE protect
    ADD COLUMN treatmentApproach varchar(500) after sourceTreatmentApproach;

#Patient db

ALTER TABLE sample
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE snpcheck
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE metric
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE flagstat
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE somaticVariant
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE germlineVariant
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE purity
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE purityRange
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE copyNumber
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE geneCopyNumber
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE driverCatalog
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE germlineDeletion
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE structuralVariant
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svAnnotation
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svCluster
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svLink
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svDriver
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svBreakend
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE svFusion
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE structuralVariantGermline
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE signature
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE chord
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE peachCalls
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE peachGenotype
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE virusBreakend
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE virusAnnotation
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE cuppa
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE protect
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE hlaType
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;

ALTER TABLE hlaTypeDetails
    ADD COLUMN isolationBarcode varchar(255) NOT NULL after sampleId;