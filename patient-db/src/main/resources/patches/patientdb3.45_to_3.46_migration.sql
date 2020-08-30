ALTER TABLE somaticVariant
    DROP COLUMN cosmicId,
    DROP COLUMN dbsnpId,
    ADD COLUMN reported BOOLEAN NOT NULL;

ALTER TABLE structuralVariant
    CHANGE COLUMN startLinkedBy startLinkedBy varchar(1024),
    CHANGE COLUMN endLinkedBy endLinkedBy varchar(1024);