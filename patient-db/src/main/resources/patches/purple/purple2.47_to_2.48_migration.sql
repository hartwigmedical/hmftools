ALTER TABLE purity
    ADD COLUMN svTmb INT not null DEFAULT 0,
    ADD COLUMN deletedGenes INT not null DEFAULT 0,
    ADD COLUMN copyNumberSegments INT not null DEFAULT 0,
    ADD COLUMN unsupportedCopyNumberSegments INT not null DEFAULT 0,
    ADD COLUMN contamination DOUBLE PRECISION not null DEFAULT 0,
    ADD COLUMN germlineAberration varchar(255) not null DEFAULT "NONE",
    ADD COLUMN amberGender varchar(255) not null,
    CHANGE COLUMN status fitMethod varchar(255) NOT NULL;

UPDATE purity set germlineAberration = "KLINEFELTER" where gender = "MALE_KLINEFELTER";
UPDATE purity set gender = 'MALE' where gender = 'MALE_KLINEFELTER';
UPDATE purity set amberGender = gender;