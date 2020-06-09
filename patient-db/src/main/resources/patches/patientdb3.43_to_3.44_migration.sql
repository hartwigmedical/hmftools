ALTER TABLE somaticVariant
    RENAME COLUMN ploidy TO variantCopyNumber,
    RENAME COLUMN minorAllelePloidy TO minorAlleleCopyNumber;

ALTER TABLE structuralVariant
    RENAME COLUMN ploidy to junctionCopyNumber;

ALTER TABLE copyNumber
    RENAME COLUMN minorAllelePloidy TO minorAlleleCopyNumber,
    RENAME COLUMN majorAllelePloidy TO majorAlleleCopyNumber;

ALTER TABLE copyNumberGermline
    RENAME COLUMN minorAllelePloidy TO minorAlleleCopyNumber,
    RENAME COLUMN majorAllelePloidy TO majorAlleleCopyNumber;

ALTER TABLE geneCopyNumber
    RENAME COLUMN minMinorAllelePloidy to minMinorAlleleCopyNumber;