ALTER TABLE structuralVariant
ADD insertSequenceRepeatClass varchar(64),
ADD insertSequenceRepeatType varchar(64),
ADD insertSequenceRepeatOrientation tinyint,
ADD insertSequenceRepeatCoverage DOUBLE PRECISION;

ALTER TABLE structuralVariant
ADD recoveryMethod varchar(64) AFTER recovered,
ADD recoveryFilter varchar(255) AFTER recoveryMethod;
