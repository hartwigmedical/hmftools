ALTER TABLE somaticVariant
    CHANGE COLUMN effect worstEffect varchar(255) NOT NULL,
    CHANGE COLUMN codingEffect worstCodingEffect varchar(255) NOT NULL,
    ADD COLUMN genesEffected int not null AFTER gene,
    ADD COLUMN worstEffectTranscript varchar(255) NOT NULL AFTER dbsnpId,
    ADD COLUMN canonicalEffect varchar(255) NOT NULL AFTER worstCodingEffect,
    ADD COLUMN canonicalCodingEffect varchar(255) NOT NULL AFTER canonicalEffect;