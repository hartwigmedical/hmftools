ALTER TABLE somaticVariant
    ADD COLUMN mappability DOUBLE PRECISION NOT NULL AFTER loh;


ALTER TABLE somaticVariant
    ADD COLUMN type varchar(255) NOT NULL AFTER filter,
    ADD INDEX (type);