ALTER TABLE somaticVariant
    ADD recovered BOOLEAN NOT NULL DEFAULT 0 AFTER minorAllelePloidy;
