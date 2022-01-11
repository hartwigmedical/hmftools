####
# Sage

ALTER TABLE somaticVariant
	DROP COLUMN phasedInframeIndel;

####
# virus-interpreter
ALTER TABLE virusAnnotation
    ADD COLUMN lilikelihood VARCHAR(10) NOT NULL,;