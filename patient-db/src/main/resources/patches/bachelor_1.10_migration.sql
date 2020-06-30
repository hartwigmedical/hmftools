ALTER TABLE germlineVariant
    CHANGE minorAllelePloidy minorAlleleJunctionCopyNumber DOUBLE PRECISION,
    DROP COLUMN cosmicId,
    DROP COLUMN dbsnpId,
    ADD COLUMN clinvarInfo VARCHAR(255);
