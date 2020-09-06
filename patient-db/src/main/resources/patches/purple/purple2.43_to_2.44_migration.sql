ALTER TABLE purity
    MODIFY COLUMN tml INT not null;

ALTER TABLE copyNumber
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION not null,
    CHANGE COLUMN majorAllelePloidy majorAlleleCopyNumber DOUBLE PRECISION not null;

ALTER TABLE copyNumberGermline
    CHANGE COLUMN minorAllelePloidy minorAlleleCopyNumber DOUBLE PRECISION not null,
    CHANGE COLUMN majorAllelePloidy majorAlleleCopyNumber DOUBLE PRECISION not null;

ALTER TABLE geneCopyNumber
    CHANGE COLUMN minMinorAllelePloidy minMinorAlleleCopyNumber DOUBLE PRECISION not null;