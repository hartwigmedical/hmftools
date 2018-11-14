ALTER TABLE structuralVariantBreakend
    ADD sampleId varchar(255) NOT NULL AFTER modified;

<<<<<<< HEAD
ALTER TABLE structuralVariantDisruption
    ADD modified DATETIME NOT NULL DEFAULT '2018-01-01' AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantFusion
    ADD modified DATETIME NOT NULL DEFAULT '2018-01-01' AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;
=======
ALTER TABLE structuralVariantFusion
    ADD modified DATETIME AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;

ALTER TABLE structuralVariantDisruption
    ADD modified DATETIME AFTER id,
    ADD sampleId varchar(255) NOT NULL AFTER modified;
>>>>>>> 96fb392a1136c412b96107186d1735be6e5ecb55
