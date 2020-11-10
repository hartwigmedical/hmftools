ALTER TABLE amberAnonymous
    ADD COLUMN deleted BOOLEAN NOT NULL;

ALTER TABLE doidNode
    ADD COLUMN snomedId varchar(255) AFTER doidTerm;