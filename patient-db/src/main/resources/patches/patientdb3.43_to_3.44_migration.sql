ALTER TABLE somaticVariant
    CHANGE COLUMN ploidy variantCopyNumber DOUBLE PRECISION NOT NULL,
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION NOT NULL;

ALTER TABLE structuralVariant
    CHANGE COLUMN ploidy junctionCopyNumber DOUBLE PRECISION;

ALTER TABLE copyNumber
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION not null,
    CHANGE COLUMN majorAllelePloidy majorAlleleCopyNumber DOUBLE PRECISION not null;

ALTER TABLE copyNumberGermline
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION not null,
    CHANGE COLUMN majorAllelePloidy majorAlleleCopyNumber DOUBLE PRECISION not null;

ALTER TABLE geneCopyNumber
    CHANGE COLUMN minMinorAllelePloidy minMinorAlleleCopyNumber DOUBLE PRECISION not null;