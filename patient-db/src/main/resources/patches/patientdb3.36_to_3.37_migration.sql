ALTER TABLE patient
    ADD blacklisted BOOLEAN NOT NULL DEFAULT FALSE after patientIdentifier;
