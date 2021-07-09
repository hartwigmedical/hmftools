DROP TABLE rna;

ALTER TABLE canonicalTranscript
    DROP COLUMN transcriptVersion;

ALTER TABLE geneCopyNumber
    DROP COLUMN transcriptVersion;
