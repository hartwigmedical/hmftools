####
# SQL updates for Pipeline release 5.28
# NOTE: only add updates to this script if the tools impacted by them will be released with this pipeline release

# SAGE
ALTER TABLE somaticVariant
	DROP COLUMN localRealignmentSet,
	CHANGE COLUMN type type varchar(10) NOT NULL,
	CHANGE COLUMN gene gene varchar(50) NOT NULL,
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50),
	CHANGE COLUMN canonicalHgvsCodingImpact canonicalHgvsCodingImpact varchar(100),
	CHANGE COLUMN canonicalHgvsProteinImpact canonicalHgvsProteinImpact varchar(200),
	CHANGE COLUMN germlineStatus germlineStatus varchar(50) NOT NULL;

ALTER TABLE germlineVariant
	CHANGE COLUMN type type varchar(10) NOT NULL,
	CHANGE COLUMN gene gene varchar(50) NOT NULL,
	CHANGE COLUMN localPhaseSet localPhaseSet varchar(50),
	CHANGE COLUMN canonicalCodingEffect canonicalCodingEffect varchar(50) NOT NULL,
    CHANGE COLUMN canonicalHgvsCodingImpact canonicalHgvsCodingImpact varchar(100),
    CHANGE COLUMN canonicalHgvsProteinImpact canonicalHgvsProteinImpact varchar(200),
    CHANGE COLUMN worstCodingEffect worstCodingEffect varchar(50);

ALTER TABLE purity
	CHANGE COLUMN version version VARCHAR(10) NOT NULL,
    CHANGE COLUMN gender gender varchar(20) NOT NULL,
    CHANGE COLUMN fitMethod fitMethod varchar(20) NOT NULL,
    ADD COLUMN runMode VARCHAR(20) NOT NULL after fitMethod,
    CHANGE COLUMN amberGender amberGender varchar(20) NOT NULL,
    ADD COLUMN targeted BOOLEAN not null;
