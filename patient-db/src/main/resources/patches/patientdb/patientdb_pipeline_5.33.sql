####
# SQL updates for Pipeline release 5.33
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Tumor-only fields

ALTER TABLE somaticVariant
    ADD COLUMN clinvarInfo varchar(255) NULL after reported,
    ADD COLUMN gnomadFrequency DOUBLE PRECISION NULL after clinvarInfo,
    ADD COLUMN somaticLikelihood varchar(10) NULL after gnomadFrequency;

ALTER TABLE purity
    DROP COLUMN version;

ALTER TABLE sample
    ADD COLUMN cohortId varchar(255) NOT NULL after sampleId
