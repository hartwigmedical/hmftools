UPDATE chord SET predictedResponse = hrd WHERE hrd > 0;

ALTER TABLE chord MODIFY predictedResponse BOOLEAN NOT NULL;

