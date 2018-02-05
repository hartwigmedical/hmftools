ALTER TABLE somaticVariant
    ADD COLUMN codingEffect varchar(255) NOT NULL AFTER effect,
    CHANGE COLUMN loh biallelic BOOLEAN NOT NULL,
    ADD COLUMN hotspot BOOLEAN NOT NULL AFTER biallelic;


ALTER TABLE geneCopyNumber
    ADD COLUMN nonsenseBiallelicVariants int not null ,
    ADD COLUMN nonsenseNonBiallelicVariants int not null,
    ADD COLUMN nonsenseNonBiallelicPloidy DOUBLE PRECISION not null,
    ADD COLUMN spliceBiallelicVariants int not null,
    ADD COLUMN spliceNonBiallelicVariants int not null,
    ADD COLUMN spliceNonBiallelicPloidy DOUBLE PRECISION not null,
    ADD COLUMN missenseBiallelicVariants int not null,
    ADD COLUMN missenseNonBiallelicVariants int not null,
    ADD COLUMN missenseNonBiallelicPloidy DOUBLE PRECISION not null,
    ADD COLUMN minMinorAllelePloidy DOUBLE PRECISION not null,
    ADD COLUMN exonicBases int not null;