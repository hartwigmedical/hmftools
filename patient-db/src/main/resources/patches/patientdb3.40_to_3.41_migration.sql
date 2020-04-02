ALTER TABLE somaticVariant
    ADD tier varchar(20) NOT NULL AFTER highConfidence;

UPDATE somaticVariant
    set tier = 'UNKNOWN';
