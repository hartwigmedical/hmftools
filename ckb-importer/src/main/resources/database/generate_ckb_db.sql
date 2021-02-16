SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS ckbEntry;

DROP TABLE IF EXISTS clinicalTrialContact;
DROP TABLE IF EXISTS clinicalTrialLocation;
DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetailMolecularProfile;
DROP TABLE IF EXISTS clinicalTrialVariantRequirementDetail;
DROP TABLE IF EXISTS clinicalTrialIndication;
DROP TABLE IF EXISTS clinicalTrialTherapy;
DROP TABLE IF EXISTS clinicalTrialAgeGroup;
DROP TABLE IF EXISTS clinicalTrial;

DROP TABLE IF EXISTS drugGlobalApproavalStatusMolecularProfile;
DROP TABLE IF EXISTS drugGlobalApproavalStatusIndication;
DROP TABLE IF EXISTS drugGlobalApproavalStatusTherapy;
DROP TABLE IF EXISTS drugGlobalApproavalStatus;
DROP TABLE IF EXISTS drugTherapy;
DROP TABLE IF EXISTS drugEvidenceTreatmentApproch;
DROP TABLE IF EXISTS drugEvidenceReference;
DROP TABLE IF EXISTS drugEvidenceIndication;
DROP TABLE IF EXISTS drugEvidenceTherapy;
DROP TABLE IF EXISTS drugEvidenceMolecularProfile;
DROP TABLE IF EXISTS drugEvidence;
DROP TABLE IF EXISTS drugClinicalTrialTherapy;
DROP TABLE IF EXISTS drugClinicalTrial;
DROP TABLE IF EXISTS drugDrugClass;
DROP TABLE IF EXISTS drugDescriptionReference;
DROP TABLE IF EXISTS drugDescription;
DROP TABLE IF EXISTS drugSynonym;
DROP TABLE IF EXISTS drugTerm;
DROP TABLE IF EXISTS drug;

SET FOREIGN_KEY_CHECKS = 1;