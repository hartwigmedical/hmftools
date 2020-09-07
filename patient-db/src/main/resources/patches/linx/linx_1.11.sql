ALTER TABLE svCluster
    CHANGE resolvedType category VARCHAR(20),
    CHANGE subType resolvedType VARCHAR(20),
    DROP COLUMN subClonal;

