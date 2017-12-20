ALTER TABLE somaticVariant
    ADD COLUMN mappability DOUBLE PRECISION NOT NULL AFTER loh,
    ADD COLUMN germlineStatus varchar(255) NOT NULL AFTER mappability,
    ADD COLUMN type varchar(255) NOT NULL AFTER filter,
    ADD INDEX (type);