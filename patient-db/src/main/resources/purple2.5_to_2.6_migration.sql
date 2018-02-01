ALTER TABLE somaticVariant
    ADD COLUMN codingEffect varchar(255) NOT NULL AFTER effect,
    CHANGE COLUMN loh biallelic BOOLEAN NOT NULL,
    ADD COLUMN hotspot BOOLEAN NOT NULL AFTER biallelic;