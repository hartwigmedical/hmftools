ALTER TABLE somaticVariant
    CHANGE adjustedCopyNumber copyNumber DOUBLE PRECISION NOT NULL,
    ADD ploidy DOUBLE PRECISION NOT NULL AFTER adjustedVaf;

UPDATE somaticVariant
SET ploidy = copyNumber * adjustedVaf;