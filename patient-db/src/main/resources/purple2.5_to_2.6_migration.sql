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
    ADD COLUMN minMinorAllelePloidy DOUBLE PRECISION not null;

ALTER TABLE geneCopyNumber
        DROP COLUMN meanCopyNumber;

DROP TABLE IF EXISTS canonicalTranscript;
CREATE TABLE canonicalTranscript
(   id int NOT NULL AUTO_INCREMENT,
    modified DATETIME NOT NULL,
    gene varchar(255) NOT NULL,
    geneId varchar(255) NOT NULL,
    chromosomeBand varchar(255) NOT NULL,
    chromosome varchar(255) NOT NULL,
    geneStart int not null,
    geneEnd int not null,
    transcriptId varchar(255) NOT NULL,
    transcriptVersion int not null,
    transcriptStart int not null,
    transcriptEnd int not null,
    exons int not null,
    exonStart int not null,
    exonEnd int not null,
    exonBases int not null,
    codingStart int not null,
    codingEnd int not null,
    codingBases int not null,
    PRIMARY KEY (id),
    INDEX(gene),
    INDEX(transcriptId)
);