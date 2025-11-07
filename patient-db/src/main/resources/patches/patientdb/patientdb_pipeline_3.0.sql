####
# SQL updates for Pipeline release 3.0
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Drivers
ALTER TABLE driverCatalog
    ADD COLUMN reportedStatus varchar(50) NULL after likelihoodMethod;
