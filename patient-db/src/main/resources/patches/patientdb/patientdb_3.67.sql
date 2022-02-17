####
# PROTECT

ALTER TABLE protect
    ADD COLUMN gene varchar(255) AFTER sampleId;

ALTER TABLE protect
    ADD COLUMN evidenceType varchar(50) NOT NULL AFTER event;

ALTER TABLE protect
    ADD COLUMN rangeRank int AFTER evidenceType;

ALTER TABLE somaticVariant
	DROP COLUMN localRealignmentSet,
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50);

ALTER TABLE germlineVariant
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50);
