####
# SQL updates for Pipeline release 3.0
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

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

