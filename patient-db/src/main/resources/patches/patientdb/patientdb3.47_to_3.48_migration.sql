ALTER TABLE baseline
    DROP COLUMN cancerSubtype,
    ADD COLUMN primaryTumorSubLocation varchar(255) AFTER primaryTumorLocation,
    ADD COLUMN primaryTumorType varchar(255) AFTER primaryTumorSubLocation,
    ADD COLUMN primaryTumorSubType varchar(255) AFTER primaryTumorType,
    ADD COLUMN primaryTumorExtraDetails varchar(255) AFTER primaryTumorSubType,
    ADD COLUMN primaryTumorOverridden BOOLEAN AFTER primaryTumorExtraDetails
