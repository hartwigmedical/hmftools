ALTER TABLE somaticVariant
    CHANGE COLUMN ploidy variantCopyNumber DOUBLE PRECISION NOT NULL,
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION NOT NULL;

ALTER TABLE structuralVariant
    CHANGE COLUMN ploidy junctionCopyNumber DOUBLE PRECISION;
