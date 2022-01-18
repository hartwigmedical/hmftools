####
# Sage

ALTER TABLE somaticVariant
	DROP COLUMN phasedInframeIndel;

####
# Virus Interpreter
ALTER TABLE virusAnnotation
    ADD COLUMN likelihood VARCHAR(10) NOT NULL AFTER reported;