####
# Sage

ALTER TABLE somaticVariant
	DROP COLUMN phasedInframeIndel;

####
# virus-interpreter
ALTER TABLE virusAnnotation
    ADD COLUMN isHighRisk BOOLEAN AFTER reported;