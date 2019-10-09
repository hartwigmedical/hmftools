ALTER TABLE chord MODIFY predictedResponse BOOLEAN NOT NULL;

update table chord set predictedResponse = hrd > 0.5;
