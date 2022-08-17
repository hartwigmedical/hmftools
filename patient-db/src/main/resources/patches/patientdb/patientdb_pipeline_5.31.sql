####
# SQL updates for Pipeline release 5.31
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# Protect

ALTER TABLE protect
    ADD COLUMN sourceTreatmentApproach varchar(500) after treatment;

ALTER TABLE protect
    ADD COLUMN treatmentApproach varchar(500) after sourceTreatmentApproach;