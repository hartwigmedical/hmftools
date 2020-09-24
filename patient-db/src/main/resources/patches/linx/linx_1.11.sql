ALTER TABLE svCluster
    CHANGE resolvedType category VARCHAR(20),
    CHANGE subType resolvedType VARCHAR(20),
    DROP COLUMN subClonal;

ALTER TABLE svFusion
    CHANGE name name VARCHAR(50) NOT NULL AFTER threePrimeBreakendId,
    CHANGE reportedType reportedType varchar(50) NOT NULL after reported,
    CHANGE phased phased VARCHAR(20) NOT NULL AFTER reportedType,
    ADD likelihood VARCHAR(10) NOT NULL AFTER phased;

