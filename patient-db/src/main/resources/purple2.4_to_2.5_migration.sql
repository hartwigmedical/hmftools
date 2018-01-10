ALTER TABLE somaticVariant
    ADD COLUMN minorAllelePloidy DOUBLE PRECISION NOT NULL AFTER germlineStatus;