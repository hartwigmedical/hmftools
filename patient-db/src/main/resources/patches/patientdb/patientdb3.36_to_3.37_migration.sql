ALTER TABLE patient
    ADD blacklisted BOOLEAN NOT NULL after patientIdentifier;
