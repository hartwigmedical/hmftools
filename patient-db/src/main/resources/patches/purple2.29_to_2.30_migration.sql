DROP TABLE IF EXISTS copyNumberRegion;

ALTER TABLE purity
  ADD COLUMN wholeGenomeDuplication BOOLEAN NOT NULL AFTER polyclonalProportion;

UPDATE purity p,
  (SELECT sampleId, COUNT(*) AS duplicatedAutosomes
      FROM
     (SELECT sampleId, chromosome, round(SUM(bafCount*baf*copyNumber)/SUM(bafCount),1) AS lwMajorAlleleAvg
    FROM copyNumber
     WHERE chromosome NOT IN ('X', 'Y')
       AND bafCount > 0
     GROUP BY sampleId, chromosome
     HAVING lwMajorAlleleAvg > 1.5) a
     GROUP BY sampleId) b
SET p.wholeGenomeDuplication = b.duplicatedAutosomes > 10
WHERE p.sampleId = b.sampleId;