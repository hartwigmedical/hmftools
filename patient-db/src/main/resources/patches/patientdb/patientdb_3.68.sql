####
# SAGE

ALTER TABLE somaticVariant
	DROP COLUMN localRealignmentSet,
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50);

ALTER TABLE germlineVariant
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50);