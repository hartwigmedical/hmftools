####
# SQL updates for Pipeline release 5.29
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Purple

ALTER TABLE germlineDeletion
    ADD COLUMN chromosomeBand VARCHAR(20) NOT NULL after chromosome;

# Protect

ALTER TABLE protect
	CHANGE COLUMN sourceUrls sourceUrls varchar(2500),
	CHANGE COLUMN evidenceUrls evidenceUrls varchar(2500);
