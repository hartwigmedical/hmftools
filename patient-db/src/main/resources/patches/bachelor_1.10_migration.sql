ALTER TABLE germlineVariant
    CHANGE minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION,
    DROP COLUMN cosmicId,
    DROP COLUMN dbsnpId,
    ADD COLUMN reported BOOLEAN NOT NULL AFTER filter,
    ADD COLUMN pathogenic VARCHAR(50) NOT NULL AFTER reported,
    ADD COLUMN clinvarInfo VARCHAR(255);
