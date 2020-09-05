ALTER TABLE canonicalTranscript
 ADD assembly varchar(255) NOT NULL AFTER modified;

UPDATE canonicalTranscript set assembly = "GRCh37";