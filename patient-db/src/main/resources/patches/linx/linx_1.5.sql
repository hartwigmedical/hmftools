ALTER TABLE svLink
    CHANGE chromosome chromosome VARCHAR(10);

ALTER TABLE structuralVariant
    DROP INDEX sampleId,
    ADD INDEX (sampleId, svId);

ALTER TABLE svBreakend
    DROP INDEX sampleId,
    ADD INDEX (sampleId, svId);

ALTER TABLE svAnnotation
    DROP INDEX sampleId,
    ADD INDEX (sampleId, svId);

ALTER TABLE svCluster
    DROP INDEX sampleId,
    ADD INDEX (sampleId, clusterId);

ALTER TABLE svLink
    DROP INDEX sampleId,
    ADD INDEX (sampleId, clusterId);

ALTER TABLE svDriver
    DROP INDEX sampleId,
    ADD INDEX (sampleId, clusterId);
