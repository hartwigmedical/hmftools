####
# SQL updates for Pipeline release 5.30
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Purple

ALTER TABLE germlineDeletion
    ADD COLUMN chromosomeBand VARCHAR(20) NOT NULL after chromosome;
