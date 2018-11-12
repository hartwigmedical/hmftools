ALTER TABLE structuralVariantBreakend
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantDisruption
    ADD modified DATETIME NOT NULL DEFAULT '2018-01-01' AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantFusion
    ADD modified DATETIME NOT NULL DEFAULT '2018-01-01' AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;
