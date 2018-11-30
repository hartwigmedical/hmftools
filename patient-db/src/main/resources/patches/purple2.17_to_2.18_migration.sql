ALTER TABLE somaticVariant
    ADD recovered BOOLEAN NOT NULL AFTER minorAllelePloidy;
