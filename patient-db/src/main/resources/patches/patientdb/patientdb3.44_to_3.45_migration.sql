ALTER TABLE somaticVariant
    DROP COLUMN highConfidence,
    ADD COLUMN localPhaseSet int,
    ADD COLUMN localRealignmentSet int,
    ADD COLUMN phasedInframeIndel int;