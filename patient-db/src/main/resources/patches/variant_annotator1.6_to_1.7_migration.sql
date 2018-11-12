ALTER TABLE structuralVariantBreakend
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantFusion
    ADD modified DATETIME AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantDisruption
    ADD modified DATETIME AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;