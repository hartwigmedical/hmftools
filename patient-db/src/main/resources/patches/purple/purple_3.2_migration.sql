ALTER TABLE geneCopyNumber
    ADD COLUMN canonicalTranscript BOOLEAN NOT NULL DEFAULT 1;

ALTER TABLE somaticVariant
    ADD COLUMN spliceRegion BOOLEAN NOT NULL DEFAULT 0 AFTER canonicalHgvsProteinImpact,
    ADD COLUMN otherTranscriptEffects varchar(255) AFTER spliceRegion,
    CHANGE COLUMN genesEffected genesAffected int not null,
    DROP COLUMN worstEffect,
    DROP COLUMN worstEffectTranscript;

ALTER TABLE driverCatalog
    ADD COLUMN transcriptId varchar(50) NOT NULL AFTER gene,
    ADD COLUMN canonicalTranscript BOOLEAN NOT NULL DEFAULT 1 after transcript;

ALTER TABLE germlineVariant
    ADD COLUMN spliceRegion BOOLEAN NOT NULL DEFAULT 0 AFTER canonicalHgvsProteinImpact,
    ADD COLUMN otherTranscriptEffects varchar(255) AFTER spliceRegion,
    CHANGE COLUMN genesEffected genesAffected int not null,
    DROP COLUMN worstEffect,
    DROP COLUMN worstEffectTranscript;
