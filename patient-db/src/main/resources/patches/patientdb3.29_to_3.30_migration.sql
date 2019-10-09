update chord set predictedResponse = hrd where hrd > 0;

ALTER TABLE chord MODIFY predictedResponse BOOLEAN NOT NULL;

