ALTER TABLE purity
    ADD COLUMN germlineAberration varchar(255) not null DEFAULT "NONE";

UPDATE purity set germlineAberration = "KLINEFELTER" where gender = "MALE_KLINEFELTER";
UPDATE purity set gender = 'MALE' where gender = 'MALE_KLINEFELTER';