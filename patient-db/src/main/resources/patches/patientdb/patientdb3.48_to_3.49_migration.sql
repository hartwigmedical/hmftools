ALTER TABLE baseline
    ADD COLUMN doid varchar(255) AFTER primaryTumorOverridden,
    ADD COLUMN doidTerm varchar(255) AFTER doid