ALTER TABLE somaticVariant
    DROP COLUMN cosmicId,
    DROP COLUMN dbsnpId,
    ADD COLUMN reported BOOLEAN NOT NULL;