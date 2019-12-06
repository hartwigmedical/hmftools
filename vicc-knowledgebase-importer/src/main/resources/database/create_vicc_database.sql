SET FOREIGN_KEY_CHECKS = 0;

DROP TABLE IF EXISTS viccEntry;
DROP TABLE IF EXISTS gene;
DROP TABLE IF EXISTS geneIdentifier;
DROP TABLE IF EXISTS featureName;
DROP TABLE IF EXISTS tag;
DROP TABLE IF EXISTS devTag;
DROP TABLE IF EXISTS feature;
DROP TABLE IF EXISTS featureInfo;
DROP TABLE IF EXISTS featureAttribute;
DROP TABLE IF EXISTS provenance;
DROP TABLE IF EXISTS synonym;
DROP TABLE IF EXISTS link;
DROP TABLE IF EXISTS sequenceOntology;
DROP TABLE IF EXISTS hierarchy;
DROP TABLE IF EXISTS association;
DROP TABLE IF EXISTS associationVariant;
DROP TABLE IF EXISTS evidence;
DROP TABLE IF EXISTS evidenceInfo;
DROP TABLE IF EXISTS evidenceType;
DROP TABLE IF EXISTS publicationUrl;
DROP TABLE IF EXISTS phenotype;
DROP TABLE IF EXISTS phenotypeType;
DROP TABLE IF EXISTS environmentalContext;
DROP TABLE IF EXISTS approvedCountry;
DROP TABLE IF EXISTS taxonomy;
DROP TABLE IF EXISTS sage;
DROP TABLE IF EXISTS oncokb;
DROP TABLE IF EXISTS oncokbBiological;
DROP TABLE IF EXISTS oncokbVariantBiological;
DROP TABLE IF EXISTS oncokbConsequencesBiological;
DROP TABLE IF EXISTS oncokbGeneBiological;
DROP TABLE IF EXISTS oncokbGeneAliasesBiological;
DROP TABLE IF EXISTS oncokbClinical;
DROP TABLE IF EXISTS oncokbDrugAbstractsClinical;
DROP TABLE IF EXISTS oncokbVariantClinical;
DROP TABLE IF EXISTS oncokbConsequencesClinical;
DROP TABLE IF EXISTS oncokbGeneClinical;
DROP TABLE IF EXISTS oncokbGeneAliasesClinical;
DROP TABLE IF EXISTS civic;
DROP TABLE IF EXISTS civicAssertions;
DROP TABLE IF EXISTS civicHGVSExpressions;
DROP TABLE IF EXISTS civicClinvarEntries;
DROP TABLE IF EXISTS civicVariantAliases;
DROP TABLE IF EXISTS civicVariantTypes;
DROP TABLE IF EXISTS civicDescription;
DROP TABLE IF EXISTS civicCoordinates;
DROP TABLE IF EXISTS civicVariantsGroups;
DROP TABLE IF EXISTS civicVariantsGroupsVariants;
DROP TABLE IF EXISTS civicVariantsGroupsCoordinates;
DROP TABLE IF EXISTS civicVariantsGroupsTypes;
DROP TABLE IF EXISTS civicEvidenceItems;
DROP TABLE IF EXISTS civicDrugs;
DROP TABLE IF EXISTS civicDisease;
DROP TABLE IF EXISTS civicEvidenceItemsSource;
DROP TABLE IF EXISTS civicEvidenceItemsPublication;
DROP TABLE IF EXISTS civicEvidenceItemsClinicalTrial;
DROP TABLE IF EXISTS civicSource;
DROP TABLE IF EXISTS civicError;
DROP TABLE IF EXISTS civicPublication;
DROP TABLE IF EXISTS civicClinicalTrial;
DROP TABLE IF EXISTS civicLifecycleActions;
DROP TABLE IF EXISTS civicLastCommentedOn;
DROP TABLE IF EXISTS civicLastCommentedOnUser;
DROP TABLE IF EXISTS civicLastCommentedOnAvatars;
DROP TABLE IF EXISTS civicLastCommentedOnOrganization;
DROP TABLE IF EXISTS civicLastCommentedOnProfileImage;
DROP TABLE IF EXISTS civicLastModified;
DROP TABLE IF EXISTS civicLastModifiedUser;
DROP TABLE IF EXISTS civicLastModifiedAvatars;
DROP TABLE IF EXISTS civicLastModifiedOrganization;
DROP TABLE IF EXISTS civicLastModifiedProfileImage;
DROP TABLE IF EXISTS civicLastReviewed;
DROP TABLE IF EXISTS civicLastReviewedUser;
DROP TABLE IF EXISTS civicLastReviewedAvatars;
DROP TABLE IF EXISTS civicLastReviewedOrganization;
DROP TABLE IF EXISTS civicLastReviewedProfileImage;
DROP TABLE IF EXISTS molecularmatch;
DROP TABLE IF EXISTS molecularmatchAst;
DROP TABLE IF EXISTS molecularmatchAstLeft;
DROP TABLE IF EXISTS molecularmatchAstRight;
DROP TABLE IF EXISTS molecularmatchAstRightRight;
DROP TABLE IF EXISTS molecularmatchAstRightLeft;
DROP TABLE IF EXISTS molecularmatchAstRightLeftRight;
DROP TABLE IF EXISTS molecularmatchAstRightLeftLeft;
DROP TABLE IF EXISTS molecularmatchInstutition;
DROP TABLE IF EXISTS molecularmatchIncludeGene0;
DROP TABLE IF EXISTS molecularmatchExternalId;
DROP TABLE IF EXISTS molecularmatchIncludeStage0;
DROP TABLE IF EXISTS molecularmatchIncludeDrug1;
DROP TABLE IF EXISTS molecularmatchIncludeCondition1;
DROP TABLE IF EXISTS molecularmatchIncludeMutation1;
DROP TABLE IF EXISTS molecularmatchIncludeCondition0;
DROP TABLE IF EXISTS molecularmatchIncludeMutation0;
DROP TABLE IF EXISTS molecularmatchCriteriaMet;
DROP TABLE IF EXISTS molecularmatchSource;
DROP TABLE IF EXISTS molecularmatchTierExplanation;
DROP TABLE IF EXISTS molecularMatchTherapeuticContext;
DROP TABLE IF EXISTS molecularMatchTags;
DROP TABLE IF EXISTS molecularMatchVariantInfo;
DROP TABLE IF EXISTS molecularMatchVariantInfoConsequences;
DROP TABLE IF EXISTS molecularMatchVariantInfoFusions;
DROP TABLE IF EXISTS molecularMatchVariantInfoLocations;
DROP TABLE IF EXISTS molecularMatchVariantInfoLocationsExonNumber;
DROP TABLE IF EXISTS molecularMatchClassification;
DROP TABLE IF EXISTS molecularMatchClassificationEnd;
DROP TABLE IF EXISTS molecularMatchClassificationStart;
DROP TABLE IF EXISTS molecularMatchClassificationChr;
DROP TABLE IF EXISTS molecularMatchClassificationPathology;
DROP TABLE IF EXISTS molecularMatchClassificationRef;
DROP TABLE IF EXISTS molecularMatchClassificationExon;
DROP TABLE IF EXISTS molecularMatchClassificationAlt;
DROP TABLE IF EXISTS molecularMatchClassificationSources;
DROP TABLE IF EXISTS molecularMatchClassificationTranscripts;
DROP TABLE IF EXISTS molecularMatchClassificationCOSMICID;
DROP TABLE IF EXISTS molecularMatchClassificationdbSNP;
DROP TABLE IF EXISTS molecularMatchClassificationPopFreqMax;
DROP TABLE IF EXISTS molecularMatchClassificationExonicFunc;
DROP TABLE IF EXISTS molecularMatchClassificationNucleotideChange;
DROP TABLE IF EXISTS molecularMatchClassificationParents;
DROP TABLE IF EXISTS molecularMatchClassificationParentsTranscripts;
DROP TABLE IF EXISTS molecularMatchCriteriaUnmet;
DROP TABLE IF EXISTS molecularMatchPrefelance;
DROP TABLE IF EXISTS molecularMatchMutations;
DROP TABLE IF EXISTS molecularMatchMutationsMutationType;
DROP TABLE IF EXISTS molecularMatchMutationsSource;
DROP TABLE IF EXISTS molecularMatchMutationsSynonyms;
DROP TABLE IF EXISTS molecularMatchMutationsPathology;
DROP TABLE IF EXISTS molecularMatchMutationscDNA;
DROP TABLE IF EXISTS molecularMatchMutationsTranscriptConsequences;
DROP TABLE IF EXISTS molecularMatchMutationsTranscriptConsequencesExonNumber;
DROP TABLE IF EXISTS molecularMatchMutationsParents;
DROP TABLE IF EXISTS molecularMatchMutationsParentsTranscript;
DROP TABLE IF EXISTS molecularMatchMutationsWGSData;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataClinvarDIS;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataClinvarSIG;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataClinvarStatus;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataClinvarDBID;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataFullAA;
DROP TABLE IF EXISTS molecularMatchMutationsWGSDataGene;
DROP TABLE IF EXISTS molecularMatchMutationsWGSMap;
DROP TABLE IF EXISTS molecularMatchMutationsWGSMapSynonyms;
DROP TABLE IF EXISTS molecularMatchMutationsWGSMapProtCoords;
DROP TABLE IF EXISTS molecularMatchMutationsGRCH37Location;
DROP TABLE IF EXISTS molecularMatchMutationsGRCH37LocationConsequences;
DROP TABLE IF EXISTS molecularMatchMutationsGRCH37LocConsequencesTxsites;
DROP TABLE IF EXISTS molecularMatchMutationsGRCH37LocConsequencesExonNumber;
DROP TABLE IF EXISTS molecularMatchMutationsFusion;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBchr;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAgene;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBtx;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAchr;
DROP TABLE IF EXISTS molecularMatchMutationsFusionIns;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBgene;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAcoord;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBori;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAband;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBband;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAori;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAtx;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBcoord;
DROP TABLE IF EXISTS molecularMatchMutationsFusionBreg;
DROP TABLE IF EXISTS molecularMatchMutationsFusionAreg;
DROP TABLE IF EXISTS molecularMatchMutationsExonsInfo;
DROP TABLE IF EXISTS molecularMatchMutationsExonsInfoBoundries;
DROP TABLE IF EXISTS molecularMatchMutationsExonsInfoBoundriesExonPosities;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon1;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon2;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon3;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon4;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon5;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon6;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon7;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon8;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon9;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon10;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon11;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon12;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon13;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon14;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon15;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon16;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon17;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon18;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon19;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon20;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon21;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon22;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon23;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon24;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon25;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon26;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon27;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon28;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon29;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon30;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon31;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon32;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon33;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon34;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon35;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon36;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon37;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon38;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon39;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon40;
DROP TABLE IF EXISTS molecularMatchMutationsExonPositiesExon41;
DROP TABLE IF EXISTS pmkb;
DROP TABLE IF EXISTS pmkbTissue;
DROP TABLE IF EXISTS pmkbTumor;
DROP TABLE IF EXISTS pmkbVariant;
DROP TABLE IF EXISTS pmkbGene;
DROP TABLE IF EXISTS molecularMatchTrials;
DROP TABLE IF EXISTS molecularMatchTrialsAlterations;
DROP TABLE IF EXISTS molecularMatchTrialsIntervations;
DROP TABLE IF EXISTS molecularMatchTrialsOtherName;
DROP TABLE IF EXISTS molecularMatchTrialsOtherGroupLabel;
DROP TABLE IF EXISTS molecularMatchTrialsOverallContact;
DROP TABLE IF EXISTS molecularMatchTrialsTags;
DROP TABLE IF EXISTS molecularMatchTrialsLocations;
DROP TABLE IF EXISTS molecularMatchTrialsContact;
DROP TABLE IF EXISTS molecularMatchTrialsGeo;
DROP TABLE IF EXISTS molecularMatchTrialsLocation;
DROP TABLE IF EXISTS molecularMatchTrialsCoordinates;
DROP TABLE IF EXISTS jaxTrials;
DROP TABLE IF EXISTS jaxTrialsIndications;
DROP TABLE IF EXISTS jaxTrialsVariantRequirementDetails;
DROP TABLE IF EXISTS jaxTrialsMolecularProfile;
DROP TABLE IF EXISTS jaxTrialsTherapies;
DROP TABLE IF EXISTS jax;
DROP TABLE IF EXISTS jaxMolecularProfile;
DROP TABLE IF EXISTS jaxTherapy;
DROP TABLE IF EXISTS jaxIndication;
DROP TABLE IF EXISTS jaxReference;
DROP TABLE IF EXISTS cgi;
DROP TABLE IF EXISTS cgicDNA;
DROP TABLE IF EXISTS cgiIndividualMutation;
DROP TABLE IF EXISTS cgigDNA;
DROP TABLE IF EXISTS cgiTranscript;
DROP TABLE IF EXISTS cgiStrand;
DROP TABLE IF EXISTS cgiInfo;
DROP TABLE IF EXISTS cgiRegion;
DROP TABLE IF EXISTS brcaPart1;
DROP TABLE IF EXISTS brcaPart2;

SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE viccEntry
(   id int NOT NULL AUTO_INCREMENT,
    source varchar(50) NOT NULL,
    PRIMARY KEY (id)
);

CREATE TABLE gene
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    geneName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE geneIdentifier
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    symbol varchar(20) NOT NULL,
    entrezId varchar(20) NOT NULL,
    ensemblGeneId varchar(20),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE featureName
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    nameOfFeature varchar(2000),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE tag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    tagName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE devTag
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    devTagName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE feature
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    name varchar(1000) NOT NULL,
    biomarkerType varchar(100),
    referenceName varchar(20),
    chromosome varchar(20),
    start varchar(20),
    end varchar(20),
    ref varchar(1000),
    alt varchar(100),
    provenanceRule varchar(20),
    geneSymbol varchar(20),
    entrezId varchar(20),
    description varchar(1000),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE featureInfo
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    germlineOrSomatic varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE featureAttribute
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    aminoAcidChange varchar(20),
    germline varchar(20),
    partnerGene varchar(20),
    description varchar(20),
    exons varchar(50),
    notes varchar(20),
    cosmic varchar(20),
    effect varchar(20),
    cnvType varchar(20),
    featureAttributeId varchar(20),
    cytoband varchar(20),
    variantType varchar(20),
    dnaChange varchar(20),
    codons varchar(50),
    chromosomeBasedCnv varchar(20),
    transcript varchar(20),
    descriptionType varchar(20),
    chromosome varchar(20),
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE provenance
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    provenanceName varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE synonym
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    synonymName varchar(150) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE link
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    linkName varchar(200) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE sequenceOntology
(   id int NOT NULL AUTO_INCREMENT,
    featureId int NOT NULL,
    soid varchar(20) NOT NULL,
    parentSoid varchar(20) NOT NULL,
    name varchar(50) NOT NULL,
    parentName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (featureId) REFERENCES feature(id)
);

CREATE TABLE hierarchy
(   id int NOT NULL AUTO_INCREMENT,
    sequenceOntologyId int NOT NULL,
    hierarchyName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (sequenceOntologyId) REFERENCES sequenceOntology(id)
);

CREATE TABLE association
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    evidenceLevel varchar(20),
    evidenceLabel varchar(20),
    responseType varchar(50),
    drugLabels varchar(1500),
    sourceLink varchar(100),
    description varchar(2500) NOT NULL,
    oncogenic varchar(100),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE associationVariant
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    variantName varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE evidence
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    description varchar(1000),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE evidenceInfo
(   id int NOT NULL AUTO_INCREMENT,
    evidenceId int NOT NULL,
    publication varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id)
);

CREATE TABLE evidenceType
(   id int NOT NULL AUTO_INCREMENT,
    evidenceId int NOT NULL,
    sourceName varchar(50) NOT NULL,
    idEvidenceType varchar(200),
    PRIMARY KEY (id),
    FOREIGN KEY (evidenceId) REFERENCES evidence(id)
);

CREATE TABLE publicationUrl
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    urlOfPublication varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE phenotype
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    description varchar(255) NOT NULL,
    family varchar(100) NOT NULL,
    idPhenotype varchar(50),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE phenotypeType
(   id int NOT NULL AUTO_INCREMENT,
    phenotypeId int NOT NULL,
    source varchar(50) NOT NULL,
    term varchar(255) NOT NULL,
    idPhenotypeType varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (phenotypeId) REFERENCES phenotype(id)
);

CREATE TABLE environmentalContext
(   id int NOT NULL AUTO_INCREMENT,
    associationId int NOT NULL,
    term varchar(255),
    description varchar(255) NOT NULL,
    source varchar(50),
    usanStem varchar(100),
    toxicity varchar(1500),
    idEnvironmentalContext varchar(20),
    PRIMARY KEY (id),
    FOREIGN KEY (associationId) REFERENCES association(id)
);

CREATE TABLE approvedCountry
(   id int NOT NULL AUTO_INCREMENT,
    environmentalContextId int NOT NULL,
    approvedCountryName varchar(20) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextId) REFERENCES environmentalContext(id)
);

CREATE TABLE taxonomy
(   id int NOT NULL AUTO_INCREMENT,
    environmentalContextId int NOT NULL,
    kingdom varchar(20) NOT NULL,
    directParent varchar(100) NOT NULL,
    class varchar(50) NOT nULL,
    subClass varchar(50),
    superClass varchar(50) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (environmentalContextId) REFERENCES environmentalContext(id)
);

CREATE TABLE sage
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    entrezId varchar(255) NOT NULL,
    clinicalManifestation varchar(255) NOT NULL,
    publicationUrl varchar(255) NOT NULL,
    germlineOrSomatic varchar(255) NOT NULL,
    evidenceLabel varchar(255) NOT NULL,
    drugLabel varchar(255) NOT NULL,
    responseType varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE oncokb
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE oncokbBiological
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    mutationEffectPmids varchar(255) NOT NULL,
    Isoform varchar(255) NOT NULL,
    entrezGeneID varchar(255) NOT NULL,
    oncogenic varchar(255) NOT NULL,
    mutationEffect varchar(255) NOT NULL,
    RefSeq varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    mutationEffectAbstracts varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE oncokbVariantBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbBiologicalId int NOT NULL,
    variantResidues varchar(255),
    proteinStart varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    proteinEnd varchar(255) NOT NULL,
    refResidues varchar(255),
    alteration varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbBiologicalId) REFERENCES oncokbBiological(id)
);

CREATE TABLE oncokbConsequencesBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbVariantBiologicalId int NOT NULL,
    term varchar(255) NOT NULL,
    description varchar(255) NOT NULL,
    isGenerallyTruncating varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbVariantBiologicalId) REFERENCES oncokbVariantBiological(id)
);

CREATE TABLE oncokbGeneBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbBiologicalId int NOT NULL,
    oncogene varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    hugoSymbol varchar(255) NOT NULL,
    curatedRefSeq varchar(255),
    entrezGeneId varchar(255) NOT NULL,
    tsg varchar(255) NOT NULL,
    curatedIsoform varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbBiologicalId) REFERENCES oncokbBiological(id)
);

CREATE TABLE oncokbGeneAliasesBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbGeneBiologicalId int NOT NULL,
    geneAliases varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbGeneBiologicalId) REFERENCES oncokbGeneBiological(id)
);

CREATE TABLE oncokbClinical
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    RefSeq varchar(255) NOT NULL,
    level varchar(255) NOT NULL,
    entrezGeneID varchar(255) NOT NULL,
    drugPmids varchar(255) NOT NULL,
    cancerType varchar(255) NOT NULL,
    drug varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    levelLabel varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE oncokbDrugAbstractsClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbClinicalId int NOT NULL,
    text varchar(255) NOT NULL,
    link varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbClinicalId) REFERENCES oncokbClinical(id)
);

CREATE TABLE oncokbVariantClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbClinicalId int NOT NULL,
    variantResidues varchar(255),
    proteinStart varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    proteinEnd varchar(255) NOT NULL,
    refResidues varchar(255),
    alteration varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbClinicalId) REFERENCES oncokbClinical(id)
);

CREATE TABLE oncokbConsequencesClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbVariantClinicalId int NOT NULL,
    term varchar(255) NOT NULL,
    description varchar(255) NOT NULL,
    isGenerallyTruncating varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbVariantClinicalId) REFERENCES oncokbVariantClinical(id)
);

CREATE TABLE oncokbGeneClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbClinicalId int NOT NULL,
    oncogene varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    hugoSymbol varchar(255) NOT NULL,
    curatedRefSeq varchar(255),
    entrezGeneId varchar(255) NOT NULL,
    tsg varchar(255) NOT NULL,
    curatedIsoform varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbClinicalId) REFERENCES oncokbClinical(id)
);

CREATE TABLE oncokbGeneAliasesClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbGeneClinicalId int NOT NULL,
    geneAliases varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbGeneClinicalId) REFERENCES oncokbGeneClinical(id)
);

CREATE TABLE civic
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    entrezName varchar(255) NOT NULL,
    civicActionabilityScore varchar(255),
    alleleRegistryId varchar(255),
    geneId varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    entrezId varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    idCivic varchar(255) NOT NULL,
    description varchar(2000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE civicAssertions
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    assertions varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicHGVSExpressions
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    hgvs_expressions varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicClinvarEntries
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    clinvarEntries varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantAliases
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    variantAliases varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantTypes
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    displayName varchar(255) NOT NULL,
    description varchar(255) NOT NULL,
    url varchar(255) NOT NULL,
    soId varchar(255) NOT NULL,
    idVariantTypes varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicDescription
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    revisionId varchar(255),
    value varchar (1000),
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicCoordinates
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    chromosome2 varchar(255),
    referenceBases varchar(255),
    start2 varchar(255),
    variantBases varchar(255),
    stop varchar(255),
    stop2 varchar(255),
    representativeTranscript2 varchar(255),
    start varchar(255),
    representativeTranscript varchar(255),
    ensemblVersion varchar(255),
    chromosome varchar(255),
    referenceBuild varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantsGroups
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    idVariantsGroups varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantsGroupsVariants
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantsGroupsId int NOT NULL,
    entrez_name varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    civic_actionability_score varchar(255),
    gene_id varchar(255) NOT NULL,
    entrez_id varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    idVariants varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantsGroupsId) REFERENCES civicVariantsGroups(id)
);

CREATE TABLE civicVariantsGroupsCoordinates
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantsGroupsVariantsId int NOT NULL,
    chromosome2 varchar(255),
    referenceBases varchar(255),
    start2 varchar(255),
    variantBases varchar(255),
    stop varchar(255),
    stop2 varchar(255),
    representativeTranscript2 varchar(255),
    start varchar(255),
    representativeTranscript varchar(255),
    ensemblVersion varchar(255),
    chromosome varchar(255),
    referenceBuild varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantsGroupsVariantsId) REFERENCES civicVariantsGroupsVariants(id)
);

CREATE TABLE civicVariantsGroupsTypes
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantsGroupsVariantsId int NOT NULL,
    displayName varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    url varchar(255) NOT NULL,
    soId varchar(255) NOT NULL,
    idVariantTypes varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantsGroupsVariantsId) REFERENCES civicVariantsGroupsVariants(id)
);

CREATE TABLE civicEvidenceItems
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    status varchar(255) NOT NULL,
    rating varchar(255),
    drugInteractionType varchar(255),
    description varchar(1500) NOT NULL,
    openChangeCount varchar(255) NOT NULL,
    evidenceType varchar(255) NOT NULL,
    variantOrigin varchar(255),
    evidenceDirection varchar(255),
    variantId varchar(255),
    clinicalSignificance varchar(255),
    evidenceLevel varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    idEvidenceItems varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicDrugs
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemsId int NOT NULL,
    pubchemId varchar(255),
    idDrugs varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemsId) REFERENCES civicEvidenceItems(id)
);

CREATE TABLE civicDisease
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemsId int NOT NULL,
    doid varchar(255),
    url varchar(255),
    displayName varchar(255) NOT NULL,
    idDisease varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemsId) REFERENCES civicEvidenceItems(id)
);

CREATE TABLE civicEvidenceItemsSource
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemsId int NOT NULL,
    status varchar(255) NOT NULL,
    openAccess varchar(255),
    name varchar(1000),
    journal varchar(255),
    citation varchar(255) NOT NULL,
    pmc_Id varchar(255),
    fullJournalTitle varchar(255),
    sourceUrl varchar(255) NOT NULL,
    pubmedId varchar(255) NOT NULL,
    isReview varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemsId) REFERENCES civicEvidenceItems(id)
);

CREATE TABLE civicEvidenceItemsPublication
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemsSourceId int NOT NULL,
    year varchar(255),
    day varchar(255),
    month varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemsSourceId) REFERENCES civicEvidenceItemsSource(id)
);

CREATE TABLE civicEvidenceItemsClinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemsSourceId int NOT NULL,
    nct_id varchar(255),
    description varchar(2000),
    clinical_trial_url varchar(255),
    name varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemsSourceId) REFERENCES civicEvidenceItemsSource(id)
);

CREATE TABLE civicSource
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    status varchar(255) NOT NULL,
    openAccess varchar(255),
    name varchar(255),
    journal varchar(255),
    citation varchar(255) NOT NULL,
    pmc_Id varchar(255),
    fullJournalTitle varchar(255),
    sourceUrl varchar(255) NOT NULL,
    pubmedId varchar(255) NOT NULL,
    isReview varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicError
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicPublication
(   id int NOT NULL AUTO_INCREMENT,
    civicSourceId int NOT NULL,
    year varchar(255),
    day varchar(255),
    month varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicSourceId) REFERENCES civicSource(id)
);

CREATE TABLE civicClinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    civicSourceId int NOT NULL,
    nct_id varchar(255),
    description varchar(1000),
    clinical_trial_url varchar(255),
    name varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicSourceId) REFERENCES civicSource(id)
);

CREATE TABLE civicLifecycleActions
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicLastCommentedOn
(   id int NOT NULL AUTO_INCREMENT,
    civicLifecycleActionsId int NOT NULL,
    timestamp varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLifecycleActionsId) REFERENCES civicLifecycleActions(id)
);

CREATE TABLE civicLastCommentedOnUser
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnId int NOT NULL,
    username varchar(255) NOT NULL,
    areaOfExpertise varchar(255) NOT NULL,
    twitterHandle varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    bio varchar(1000) NOT NULL,
    url varchar(255) NOT NULL,
    createdAt varchar(255) NOT NULL,
    acceptedLicense varchar(255),
    affiliation varchar(255) NOT NULL,
    avatarUrl varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    facebookProfile varchar(255) NOT NULL,
    linkedinProfile varchar(255) NOT NULL,
    orcid varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    lastSeenAt varchar(255) NOT NULL,
    featuredExpert varchar(255) NOT NULL,
    idUser varchar(255) NOT NULL,
    signupComplete varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnId) REFERENCES civicLastCommentedOn(id)
);

CREATE TABLE civicLastCommentedOnAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnUserId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnUserId) REFERENCES civicLastCommentedOnUser(id)
);

CREATE TABLE civicLastCommentedOnOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnUserId int NOT NULL,
    url varchar(255) NOT NULL,
    idOrganization varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnUserId) REFERENCES civicLastCommentedOnUser(id)
);

CREATE TABLE civicLastCommentedOnProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnOrganizationId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnOrganizationId) REFERENCES civicLastCommentedOnOrganization(id)
);

CREATE TABLE civicLastModified
(   id int NOT NULL AUTO_INCREMENT,
    civicLifecycleActionsId int NOT NULL,
    timestamp varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLifecycleActionsId) REFERENCES civicLifecycleActions(id)
);

CREATE TABLE civicLastModifiedUser
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedId int NOT NULL,
    username varchar(255) NOT NULL,
    areaOfExpertise varchar(255),
    twitterHandle varchar(255),
    name varchar(255) NOT NULL,
    bio varchar(1500),
    url varchar(255),
    createdAt varchar(255) NOT NULL,
    acceptedLicense varchar(255),
    affiliation varchar(255),
    avatarUrl varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    facebookProfile varchar(255),
    linkedinProfile varchar(255),
    orcid varchar(255),
    displayName varchar(255) NOT NULL,
    lastSeenAt varchar(255),
    featuredExpert varchar(255) NOT NULL,
    idUser varchar(255) NOT NULL,
    signupComplete varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedId) REFERENCES civicLastModified(id)
);

CREATE TABLE civicLastModifiedAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedId) REFERENCES civicLastModified(id)
);

CREATE TABLE civicLastModifiedOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedUserId int NOT NULL,
    url varchar(255),
    idOrganization varchar(255),
    description varchar(1000),
    name varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedUserId) REFERENCES civicLastModifiedUser(id)
);

CREATE TABLE civicLastModifiedProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedOrganizationId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedOrganizationId) REFERENCES civicLastModifiedOrganization(id)
);

CREATE TABLE civicLastReviewed
(   id int NOT NULL AUTO_INCREMENT,
    civicLifecycleActionsId int NOT NULL,
    timestamp varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLifecycleActionsId) REFERENCES civicLifecycleActions(id)
);

CREATE TABLE civicLastReviewedUser
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedId int NOT NULL,
    username varchar(255) NOT NULL,
    areaOfExpertise varchar(255) NOT NULL,
    twitterHandle varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    bio varchar(1000) NOT NULL,
    url varchar(255) NOT NULL,
    createdAt varchar(255) NOT NULL,
    acceptedLicense varchar(255),
    affiliation varchar(255) NOT NULL,
    avatarUrl varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    facebookProfile varchar(255) NOT NULL,
    linkedinProfile varchar(255) NOT NULL,
    orcid varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    lastSeenAt varchar(255) NOT NULL,
    featuredExpert varchar(255) NOT NULL,
    idUser varchar(255) NOT NULL,
    signupComplete varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedId) REFERENCES civicLastReviewed(id)
);

CREATE TABLE civicLastReviewedAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedId) REFERENCES civicLastReviewed(id)
);

CREATE TABLE civicLastReviewedOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedUserId int NOT NULL,
    url varchar(255) NOT NULL,
    idOrganization varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedUserId) REFERENCES civicLastReviewedUser(id)
);

CREATE TABLE civicLastReviewedProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedOrganizationId int NOT NULL,
    x32 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
    x14 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedOrganizationId) REFERENCES civicLastReviewedOrganization(id)
);

CREATE TABLE molecularmatch
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    score varchar(255) NOT NULL,
    autoGenerateNarrative varchar(255) NOT NULL,
    clinicalSignificance varchar(255) NOT NULL,
    idMolecularMatch varchar(255) NOT NULL,
    uniqueKey varchar(255) NOT NULL,
    civicValue varchar(255) NOT NULL,
    hashKey varchar(255) NOT NULL,
    regulatoryBodyApproved varchar(255) NOT NULL,
    version varchar(255) NOT NULL,
    guidelineBody varchar(255),
    regulatoryBody varchar(255) NOT NULL,
    customer varchar(255) NOT NULL,
    direction varchar(255) NOT NULL,
    ampcap varchar(255) NOT NULL,
    guidelineVersion varchar(255),
    tier varchar(255) NOT NULL,
    mvld varchar(255) NOT NULL,
    sixtier varchar(255) NOT NULL,
    noTherapyAvailable varchar(255),
    narrative varchar(1000) NOT NULL,
    expression varchar(1000) NOT NULL,
    biomarkerClass varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE molecularmatchAst
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchAstLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstId int NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstId) REFERENCES molecularmatchAst(id)
);

CREATE TABLE molecularmatchAstRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstId int NOT NULL,
    raw varchar(1000),
    value varchar(1000),
    operator varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstId) REFERENCES molecularmatchAst(id)
);

CREATE TABLE molecularmatchAstRightRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstRightId int NOT NULL,
    raw varchar(255),
    value varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstRightId) REFERENCES molecularmatchAstRight(id)
);

CREATE TABLE molecularmatchAstRightLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstRightId int NOT NULL,
    raw varchar(255),
    value varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstRightId) REFERENCES molecularmatchAstRight(id)
);

CREATE TABLE molecularmatchAstRightLeftRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstRightLeftId int NOT NULL,
    raw varchar(255),
    value varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstRightLeftId) REFERENCES molecularmatchAstRightLeft(id)
);

CREATE TABLE molecularmatchAstRightLeftLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchAstRightLeftId int NOT NULL,
    raw varchar(255),
    operator varchar(255),
    value varchar(255),
    type varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchAstRightLeftId) REFERENCES molecularmatchAstRightLeft(id)
);

CREATE TABLE molecularmatchInstutition
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    institution varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeGene0
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeGene0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchExternalId
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    external_id varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeStage0
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeStage0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeDrug1
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeDrug1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeCondition1
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeCondition1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeMutation1
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeMutation1 varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeCondition0
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeCondition0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchIncludeMutation0
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    includeMutation0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchCriteriaMet
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    criteriaMet varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchSource
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    name varchar(255) NOT NULL,
    suppress varchar(255) NOT NULL,
    pubId varchar(255) NOT NULL,
    subType varchar(255),
    valid varchar(255) NOT NULL,
    link varchar(255) NOT NULL,
    year varchar(255) NOT NULL,
    trialId varchar(255),
    type varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    institution varchar(255),
    trialPhase varchar(255),
    functionalConsequence varchar(255),
    trustRating varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularmatchTierExplanation
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    tier varchar(255) NOT NULL,
    step varchar(255) NOT NULL,
    message varchar(1500) NOT NULL,
    success varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchTherapeuticContext
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    facet varchar(255) NOT NULL,
    suppress varchar(255),
    valid varchar(255),
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchTags
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    priority varchar(255) NOT NULL,
    compositeKey varchar(3000),
    suppress varchar(255),
    filterType varchar(255),
    term varchar(1000) NOT NULL,
    primaryValue varchar(255),
    facet varchar(255) NOT NULL,
    valid varchar(255),
    custom varchar(255),
    isNew varchar(255),
    generatedBy varchar(1000),
    manualSuppress varchar(255),
    generatedByTerm varchar(1000),
    transcript varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchVariantInfo
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    classification varchar(255) NOT NULL,
    name varchar(1000) NOT NULL,
    geneFusionPartner varchar(255) NOT NULL,
    COSMIC_ID varchar(255),
    gene varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    popFreqMax varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchVariantInfoConsequences
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    consequences varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoFusions
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    referenceGenome varchar(255) NOT NULL,
    LBPWREP varchar(255) NOT NULL,
    RBPWREP varchar(255) NOT NULL,
    exonNumber varchar(255) NOT NULL,
    chr varchar(255) NOT NULL,
    RBPWLEP varchar(255) NOT NULL,
    intronNumber varchar(255) NOT NULL,
    LBPWLEP varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoLocations
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    aminoAcidChange varchar(1000),
    intronNumber varchar(255),
    stop varchar(255) NOT NULL,
    start varchar(255) NOT NULL,
    chr varchar(255) NOT NULL,
    strand varchar(255),
    alt varchar(255),
    referenceGenome varchar(255),
    ref varchar(1000),
    cdna varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoLocationsExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoLocationsId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoLocationsId) REFERENCES molecularMatchVariantInfoLocations(id)
);

CREATE TABLE molecularMatchClassification
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    classification varchar(255) NOT NULL,
    classificationOverride varchar(255),
    geneSymbol varchar(255),
    description varchar(1500),
    priority varchar(255),
    expandGeneSearch varchar(255),
    drugsExperimentalCount varchar(255),
    drugsApprovedOffLabelCount varchar(255),
    copyNumberType varchar(255),
    publicationCount varchar(255),
    transcript varchar(255),
    name varchar(1000) NOT NULL,
    rootTerm varchar(1000),
    drugsApprovedOnLabelCount varchar(255),
    trialCount varchar(255),
    alias varchar(1000),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchClassificationEnd
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    End varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationStart
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    Start varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationChr
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    chromosome varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationPathology
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    pathology varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationRef
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    ref varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);


CREATE TABLE molecularMatchClassificationExon
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    exon varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationAlt
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    alt varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationSources
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationTranscripts
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationCOSMICID
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    COSMIC_ID varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationdbSNP
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    dbSNP varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationPopFreqMax
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    popFreqMax varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationExonicFunc
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    exonicFunc varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationNucleotideChange
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    NucleotideChange varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationParents
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    type varchar(255),
    name varchar(255) NOT NULL,
    actionableParent varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationParentsTranscripts
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationParentsId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationParentsId) REFERENCES molecularMatchClassificationParents(id)
);

CREATE TABLE molecularMatchCriteriaUnmet
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    priority varchar(255) NOT NULL,
    compositeKey varchar(3000) NOT NULL,
    isNew varchar(255),
    generatedBy varchar(255),
    manualSuppress varchar(255),
    generatedByTerm varchar(255),
    suppress varchar(255) NOT NULL,
    filterType varchar(255) NOT NULL,
    term varchar(1000) NOT NULL,
    primaryValue varchar(255),
    facet varchar(255) NOT NULL,
    valid varchar(255),
    custom varchar(255),
    transcript varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchPrefelance
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    count varchar(255) NOT NULL,
    percent varchar(255) NOT NULL,
    studyId varchar(255) NOT NULL,
    samples varchar(255) NOT NULL,
    conditionValue varchar(255),
    molecular varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchMutations
(   id int NOT NULL AUTO_INCREMENT,
    molecularmatchId int NOT NULL,
    longestTranscript varchar(255),
    description varchar(1500) NOT NULL,
    src varchar(255) NOT NULL,
    uniprotTranscript varchar(255),
    transcriptRecognized varchar(255),
    geneSymbol varchar(255) NOT NULL,
    transcript varchar(255),
    idMutations varchar(255) NOT NULL,
    name varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularmatchId) REFERENCES molecularmatch(id)
);

CREATE TABLE molecularMatchMutationsMutationType
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    mutationType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsSource
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsSynonyms
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    synonyms varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsPathology
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    pathology varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationscDNA
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    cDNA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsTranscriptConsequences
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    aminoAcidChange varchar(255),
    compositeKey varchar(255) NOT NULL,
    intronNumber varchar(255),
    suppress varchar(255) NOT NULL,
    stop varchar(255),
    custom varchar(255) NOT NULL,
    start varchar(255),
    chr varchar(255),
    strand varchar(255) NOT NULL,
    validated varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    cdna varchar(255),
    referenceGenome varchar(255) NOT NULL,
    ref varchar(1000),
    alt varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsTranscriptConsequencesExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsTranscriptConsequencesId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsTranscriptConsequencesId) REFERENCES molecularMatchMutationsTranscriptConsequences(id)
);

CREATE TABLE molecularMatchMutationsParents
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    type varchar(255),
    name varchar(255) NOT NULL,
    actionableParent varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsParentsTranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsParentsId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsParentsId) REFERENCES molecularMatchMutationsParents(id)
);

CREATE TABLE molecularMatchMutationsWGSData
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    exonicFunc varchar(255),
    dbSNP varchar(255),
    ExAC_NFE varchar(255),
    ExAC_FIN varchar(255),
    G1000_ALL varchar(255),
    G1000_SAS varchar(255),
    G1000_EAS varchar(255),
    G1000_AFR varchar(255),
    ExAC_SAS varchar(255),
    ExAC_EAS varchar(255),
    ExAC_AMR varchar(255),
    ExAC_AFR varchar(255),
    ExAC_Freq varchar(255),
    End varchar(255) NOT NULL,
    start varchar(255) NOT NULL,
    SiPhy_29way_logOdds varchar(255) NOT NULL,
    Ref varchar(255) NOT NULL,
    GERP_RS varchar(255) NOT NULL,
    FATHMM varchar(255) NOT NULL,
    NucleotideChange varchar(255) NOT NULL,
    phyloP100way_vertebrate varchar(255) NOT NULL,
    Func varchar(255) NOT NULL,
    GWAS_PUBMED varchar(255),
    Transcript varchar(255) NOT NULL,
    ESP6500si_AA varchar(255),
    ESP6500si_EA varchar(255),
    G1000_EUR varchar(255),
    G1000_AMR varchar(255),
    Chr_Start_Ref_Alt varchar(255) NOT NULL,
    AA varchar(255),
    PopFreqMax varchar(255) NOT NULL,
    FATHMM_Pred varchar(255) NOT NULL,
    wgRna varchar(255),
    phyloP46way_placental varchar(255) NOT NULL,
    keyValue varchar(255) NOT NULL,
    targetScanS varchar(255),
    Chr varchar(255) NOT NULL,
    COSMIC_ID varchar(255),
    alt varchar(255) NOT NULL,
    GWAS_DIS varchar(255),
    GWAS_SNP varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsWGSDataClinvarDIS
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    clinvarDIS varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSDataClinvarSIG
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    clinvarSIG varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSDataClinvarStatus
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    clinvarStatus varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSDataClinvarDBID
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    clinvarDBID varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSDataFullAA
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    fullAA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSDataGene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSDataId int NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSDataId) REFERENCES molecularMatchMutationsWGSData(id)
);

CREATE TABLE molecularMatchMutationsWGSMap
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    AA varchar(255),
    name varchar(255) NOT NULL,
    GRCh37_Chr_Start_Ref_Alt varchar(255) NOT NULL,
    NucleotideChange varchar(255) NOT NULL,
    Exon varchar(255),
    Gene varchar(255) NOT NULL,
    Transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsWGSMapSynonyms
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSMapId int NOT NULL,
    synonyms varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSMapId) REFERENCES molecularMatchMutationsWGSMap(id)
);

CREATE TABLE molecularMatchMutationsWGSMapProtCoords
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsWGSMapId int NOT NULL,
    protCoords varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsWGSMapId) REFERENCES molecularMatchMutationsWGSMap(id)
);

CREATE TABLE molecularMatchMutationsGRCH37Location
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    compositeKey varchar(255) NOT NULL,
    ref varchar(1000),
    stop varchar(255),
    start varchar(255),
    chr varchar(255),
    alt varchar(255),
    validated varchar(255) NOT NULL,
    strand varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsGRCH37LocationConsequences
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsGRCH37LocationId int NOT NULL,
    aminoAcidChange varchar(255),
    intronNumber varchar(255),
    transcript varchar(255) NOT NULL,
    cdna varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsGRCH37LocationId) REFERENCES molecularMatchMutationsGRCH37Location(id)
);

CREATE TABLE molecularMatchMutationsGRCH37LocConsequencesTxsites
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsGRCH37LocationConsequencesId int NOT NULL,
    txSites varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsGRCH37LocationConsequencesId) REFERENCES molecularMatchMutationsGRCH37LocationConsequences(id)
);

CREATE TABLE molecularMatchMutationsGRCH37LocConsequencesExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsGRCH37LocationConsequencesId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsGRCH37LocationConsequencesId) REFERENCES molecularMatchMutationsGRCH37LocationConsequences(id)
);

CREATE TABLE molecularMatchMutationsFusion
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    synonym varchar(255),
    source varchar(255),
    Paper varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsFusionBchr
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    Brchr varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAgene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    Agene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBtx
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    btx varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAchr
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    achr varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionIns
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    ins varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBgene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    bgene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAcoord
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    acoord varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBori
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    bori varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAband
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    aband varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBband
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    bband varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAori
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    aori varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAtx
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    atx varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBcoord
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    bcoord varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionBreg
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    num varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsFusionAreg
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsFusionId int NOT NULL,
    num varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsFusionId) REFERENCES molecularMatchMutationsFusion(id)
);

CREATE TABLE molecularMatchMutationsExonsInfo
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsId int NOT NULL,
    txStart varchar(255),
    cdsEnd varchar(255),
    chr varchar(255) NOT NULL,
    cdsStart varchar(255),
    transcript varchar(255) NOT NULL,
    txEnd varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsId) REFERENCES molecularMatchMutations(id)
);

CREATE TABLE molecularMatchMutationsExonsInfoBoundries
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonsInfoId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonsInfoId) REFERENCES molecularMatchMutationsExonsInfo(id)
);

CREATE TABLE molecularMatchMutationsExonsInfoBoundriesExonPosities
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonsInfoBoundriesId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonsInfoBoundriesId) REFERENCES molecularMatchMutationsExonsInfoBoundries(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon2
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon3
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon4
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon5
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon6
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon7
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon8
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon9
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon10
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon11
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon12
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon13
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon14
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon15
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon16
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon17
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon18
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon19
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon20
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon21
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon22
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon23
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon24
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon25
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon26
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon27
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon28
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon29
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon30
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon31
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon32
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon33
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon34
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon35
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon36
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon37
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon38
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon39
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon40
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE molecularMatchMutationsExonPositiesExon41
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationsExonPositiesId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationsExonPositiesId) REFERENCES molecularMatchMutationsExonsInfoBoundriesExonPosities(id)
);

CREATE TABLE pmkb
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE pmkbTissue
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    idTissue varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE pmkbTumor
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    idTumor varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE pmkbVariant
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    aminoAcidChange varchar(255),
    germline varchar(255),
    partnerGene varchar(255),
    codons varchar(255),
    description varchar(255),
    exons varchar(255),
    notes varchar(255),
    cosmic varchar(255),
    effect varchar(255),
    cnvType varchar(255),
    idVariant varchar(255),
    cytoband varchar(255),
    variantType varchar(255),
    dnaChange varchar(255),
    coordinates varchar(255),
    chromosomeBasedCnv varchar(255),
    transcript varchar(255),
    descriptionType varchar(255),
    chromosome varchar(255),
    name varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE pmkbGene
(   id int NOT NULL AUTO_INCREMENT,
    pmkbVariantId int NOT NULL,
    description varchar(255),
    createdAt varchar(255) NOT NULL,
    updatedAt varchar(255) NOT NULL,
    activeInd varchar(255) NOT NULL,
    externalId varchar(255) NOT NULL,
    idGene varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (pmkbVariantId) REFERENCES pmkbVariant(id)
);

CREATE TABLE molecularMatchTrials
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    status varchar(500) NOT NULL,
    startDate varchar(255),
    title varchar(1000) NOT NULL,
    score varchar(255) NOT NULL,
    briefTitle varchar(1000),
    link varchar(255) NOT NULL,
    phase varchar(255) NOT NULL,
    idMolecularMatchTrials varchar(255) NOT NULL,
    studyType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE molecularMatchTrialsAlterations
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    molecularAlterations varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsIntervations
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    intervention_name varchar(9000),
    description varchar(2000),
    intervention_type varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsOtherName
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsIntervationsId int NOT NULL,
    other_name varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsIntervationsId) REFERENCES molecularMatchTrialsIntervations(id)
);

CREATE TABLE molecularMatchTrialsOtherGroupLabel
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsIntervationsId int NOT NULL,
    arm_group_label varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsIntervationsId) REFERENCES molecularMatchTrialsIntervations(id)
);

CREATE TABLE molecularMatchTrialsOverallContact
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    phone varchar(250),
    last_name varchar(250),
    phone_ext varchar(250),
    country varchar(250),
    email varchar(250),
    affiliation varchar(250),
    city varchar(250),
    name varchar(250),
    zip varchar(250),
    url varchar(250),
    street varchar(250),
    type varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsTags
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    facet varchar(250) NOT NULL,
    compositeKey varchar(1000) NOT NULL,
    suppress varchar(250) NOT NULL,
    filterType varchar(250) NOT NULL,
    term varchar(1000) NOT NULL,
    custom varchar(250) NOT NULL,
    priority varchar(250) NOT NULL,
    alias varchar(1000),
    manualSuppress varchar(250),
    generatedBy varchar(250),
    generatedByTerm varchar(250),
    idTags varchar(250),
    manualPriority varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsLocations
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    status varchar(250) NOT NULL,
    last_name varchar(1000),
    email varchar(1000),
    phone varchar(250),
    phone_backup varchar(250),
    email_backup varchar(250),
    last_name_backup varchar(250),
    phone_ext_backup varchar(250),
    phone_ext varchar(250),
    city varchar(250),
    valid varchar(250),
    zip varchar(250),
    created varchar(250),
    country varchar(250),
    number varchar(250),
    idLocations varchar(250),
    lastUpdated varchar(250),
    state varchar(250),
    street varchar(250),
    po_box varchar(250),
    failedGeocode varchar(250),
    validMessage varchar(250),
    name varchar(1000),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsContact
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationsId int NOT NULL,
    phone varchar(250),
    name varchar(250),
    email varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationsId) REFERENCES molecularMatchTrialsLocations(id)
);

CREATE TABLE molecularMatchTrialsGeo
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationsId int NOT NULL,
    lat varchar(250),
    lon varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationsId) REFERENCES molecularMatchTrialsLocations(id)
);

CREATE TABLE molecularMatchTrialsLocation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationsId int NOT NULL,
    type varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationsId) REFERENCES molecularMatchTrialsLocations(id)
);

CREATE TABLE molecularMatchTrialsCoordinates
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationId int NOT NULL,
    coordinates varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationId) REFERENCES molecularMatchTrialsLocation(id)
);

CREATE TABLE jaxTrials
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    title varchar(500) NOT NULL,
    gender varchar(255),
    nctId varchar(255) NOT NULL,
    sponsors varchar(255) NOT NULL,
    recruitment varchar(255) NOT NULL,
    variantRequirements varchar(255) NOT NULL,
    updateDate varchar(255) NOT NULL,
    phase varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE jaxTrialsIndications
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    source varchar(255) NOT NULL,
    idIndications varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsId) REFERENCES jaxTrials(id)
);

CREATE TABLE jaxTrialsVariantRequirementDetails
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    requirementType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsId) REFERENCES jaxTrials(id)
);

CREATE TABLE jaxTrialsMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsVariantRequirementDetailsId int NOT NULL,
    profileName varchar(255) NOT NULL,
    idMolecularProfile varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsVariantRequirementDetailsId) REFERENCES jaxTrialsVariantRequirementDetails(id)
);

CREATE TABLE jaxTrialsTherapies
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    idTherapies varchar(255) NOT NULL,
    therapyName varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsId) REFERENCES jaxTrials(id)
);

CREATE TABLE jax
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    responseType varchar(255) NOT NULL,
    approvalStatus varchar(255) NOT NULL,
    evidenceType varchar(255) NOT NULL,
    efficacyEvidence varchar(1000) NOT NULL,
    idJaxEntry varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE jaxMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    jaxId int NOT NULL,
    profileName varchar(255) NOT NULL,
    idMolecularProfile varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxId) REFERENCES jax(id)
);

CREATE TABLE jaxTherapy
(   id int NOT NULL AUTO_INCREMENT,
    jaxId int NOT NULL,
    therapyName varchar(255) NOT NULL,
    idTherapy varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxId) REFERENCES jax(id)
);

CREATE TABLE jaxIndication
(   id int NOT NULL AUTO_INCREMENT,
    jaxId int NOT NULL,
    source varchar(255) NOT NULL,
    idIndication varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxId) REFERENCES jax(id)
);

CREATE TABLE jaxReference
(   id int NOT NULL AUTO_INCREMENT,
    jaxId int NOT NULL,
    url varchar(255),
    idReference varchar(255),
    pubMedId varchar(255),
    title varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (jaxId) REFERENCES jax(id)
);

CREATE TABLE cgi
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    targeting varchar(255) NOT NULL,
    source varchar(255) NOT NULL,
    primaryTumorType varchar(255) NOT NULL,
    drugsFullName varchar(255) NOT NULL,
    curator varchar(255) NOT NULL,
    drugFamily varchar(255) NOT NULL,
    alteration varchar(1000) NOT NULL,
    drug varchar(255) NOT NULL,
    biomarker varchar(1000) NOT NULL,
    drugStatus varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    assayType varchar(255) NOT NULL,
    alterationType varchar(255) NOT NULL,
    evidenceLevel varchar(255) NOT NULL,
    association varchar(255) NOT NULL,
    metastaticTumorType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE cgicDNA
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    cDNA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiIndividualMutation
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    individualMutation varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgigDNA
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    gDNA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiTranscript
( id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiStrand
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    strand varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiInfo
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    info varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiRegion
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    region varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE brcaPart1
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    Variant_frequency_LOVD TINYTEXT NOT NULL,
    Allele_frequency_FIN_ExAC TINYTEXT NOT NULL,
    ClinVarAccession_ENIGMA TINYTEXT NOT NULL,
    Homozygous_count_AFR_ExAC TINYTEXT NOT NULL,
    BX_ID_ExAC TINYTEXT NOT NULL,
    Variant_in_LOVD TINYTEXT NOT NULL,
    Allele_frequency_AFR_ExAC TINYTEXT NOT NULL,
    DBID_LOVD TINYTEXT NOT NULL,
    Chr TINYTEXT NOT NULL,
    BX_ID_ENIGMA TINYTEXT NOT NULL,
    Co_occurrence_LR_exLOVD TINYTEXT NOT NULL,
    Homozygous_count_EAS_ExAC TINYTEXT NOT NULL,
    Submitter_ClinVar TEXT NOT NULL,
    Allele_frequency_EAS_ExAC TINYTEXT NOT NULL,
    Hg37_End TINYTEXT NOT NULL,
    Submitters_LOVD TEXT NOT NULL,
    Clinical_classification_BIC TINYTEXT NOT NULL,
    Homozygous_count_NFE_ExAC TINYTEXT NOT NULL,
    Allele_count_SAS_ExAC TINYTEXT NOT NULL,
    Method_ClinVar TINYTEXT NOT NULL,
    Allele_count_NFE_ExAC TINYTEXT NOT NULL,
    Pathogenicity_all TINYTEXT NOT NULL,
    Germline_or_Somatic_BIC TINYTEXT NOT NULL,
    Homozygous_count_SAS_ExAC TINYTEXT NOT NULL,
    BIC_Nomenclature TINYTEXT NOT NULL,
    Assertion_method_ENIGMA TINYTEXT NOT NULL,
    Literature_source_exLOVD TINYTEXT NOT NULL,
    Change_Type_id TINYTEXT NOT NULL,
    Collection_method_ENIGMA TINYTEXT NOT NULL,
    Sum_family_LR_exLOVD TINYTEXT NOT NULL,
    HGVS_cDNA_LOVD TINYTEXT NOT NULL,
    Homozygous_count_FIN_ExAC TINYTEXT NOT NULL,
    EAS_Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    Ethnicity_BIC TEXT NOT NULL,
    Individuals_LOVD TINYTEXT NOT NULL,
    Variant_in_ExAC TINYTEXT NOT NULL,
    URL_ENIGMA TINYTEXT NOT NULL,
    Allele_Origin_ClinVar TINYTEXT NOT NULL,
    Allele_frequency_AMR_ExAC TINYTEXT NOT NULL,
    Variant_in_1000_Genomes TINYTEXT NOT NULL,
    AFR_Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    BX_ID_exLOVD TINYTEXT NOT NULL,
    Source TINYTEXT NOT NULL,
    Condition_ID_value_ENIGMA TINYTEXT NOT NULL,
    HGVS_Protein TINYTEXT NOT NULL,
    Ref TINYTEXT NOT NULL,
    Allele_number_AFR_ExAC TINYTEXT NOT NULL,
    Allele_count_AFR_ExAC TINYTEXT NOT NULL,
    BX_ID_LOVD TINYTEXT NOT NULL,
    Synonyms TEXT NOT NULL,
    Gene_Symbol TINYTEXT NOT NULL,
    Comment_on_clinical_significance_ENIGMA TEXT NOT NULL,
    Missense_analysis_prior_probability_exLOVD TINYTEXT NOT NULL,
    Allele_number_FIN_ExAC TINYTEXT NOT NULL,
    Posterior_probability_exLOVD TINYTEXT NOT NULL,
    Polyphen_Score TINYTEXT NOT NULL,
    Reference_Sequence TINYTEXT NOT NULL,
    Allele_count_EAS_ExAC TINYTEXT NOT NULL,
    Hg38_End TINYTEXT NOT NULL,
    HGVS_cDNA TINYTEXT NOT NULL,
    Functional_analysis_technique_LOVD TINYTEXT NOT NULL,
    SAS_Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    RNA_LOVD TINYTEXT NOT NULL,
    Combined_prior_probablility_exLOVD TINYTEXT NOT NULL,
    BX_ID_ClinVar TINYTEXT NOT NULL,
    IARC_class_exLOVD TINYTEXT NOT NULL,
    BX_ID_BIC varchar(12500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE brcaPart2
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    Sift_Prediction TINYTEXT NOT NULL,
    Allele_number_NFE_ExAC TINYTEXT NOT NULL,
    Allele_origin_ENIGMA TINYTEXT NOT NULL,
    Allele_number_OTH_ExAC TINYTEXT NOT NULL,
    Hg36_End TINYTEXT NOT NULL,
    Allele_frequency_SAS_ExAC TINYTEXT NOT NULL,
    Date_Last_Updated_ClinVar TEXT NOT NULL,
    Allele_number_EAS_ExAC TINYTEXT NOT NULL,
    Allele_frequency_OTH_ExAC TINYTEXT NOT NULL,
    Source_URL TEXT NOT NULL,
    SCV_ClinVar TEXT NOT NULL,
    Pathogenicity_expert TINYTEXT NOT NULL,
    Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    Functional_analysis_result_LOVD TINYTEXT NOT NULL,
    AMR_Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    Variant_in_ESP TINYTEXT NOT NULL,
    Variant_in_BIC TINYTEXT NOT NULL,
    Clinical_significance_ENIGMA TINYTEXT NOT NULL,
    Max_Allele_Frequency TINYTEXT NOT NULL,
    Allele_count_AMR_ExAC TINYTEXT NOT NULL,
    Variant_in_ENIGMA TINYTEXT NOT NULL,
    BX_ID_ESP TINYTEXT NOT NULL,
    Patient_nationality_BIC TEXT NOT NULL,
    BX_ID_1000_Genomes TINYTEXT NOT NULL,
    Genomic_Coordinate_hg37 TINYTEXT NOT NULL,
    Genomic_Coordinate_hg36 TINYTEXT NOT NULL,
    EUR_Allele_frequency_1000_Genomes TINYTEXT NOT NULL,
    Number_of_family_member_carrying_mutation_BIC TINYTEXT NOT NULL,
    Segregation_LR_exLOVD TINYTEXT NOT NULL,
    Allele_Frequency TINYTEXT NOT NULL,
    Minor_allele_frequency_percent_ESP TINYTEXT NOT NULL,
    Allele_frequency_ExAC TINYTEXT NOT NULL,
    Mutation_type_BIC TINYTEXT NOT NULL,
    Assertion_method_citation_ENIGMA TINYTEXT NOT NULL,
    Condition_ID_type_ENIGMA TINYTEXT NOT NULL,
    Allele_count_OTH_ExAC TINYTEXT NOT NULL,
    HGVS_protein_LOVD TINYTEXT NOT NULL,
    Variant_in_ClinVar TINYTEXT NOT NULL,
    Clinical_importance_BIC TINYTEXT NOT NULL,
    Discordant TINYTEXT NOT NULL,
    Allele_count_FIN_ExAC TINYTEXT NOT NULL,
    Condition_category_ENIGMA TINYTEXT NOT NULL,
    Allele_Frequency_ESP TINYTEXT NOT NULL,
    Homozygous_count_OTH_ExAC TINYTEXT NOT NULL,
    Genetic_origin_LOVD TINYTEXT NOT NULL,
    Homozygous_count_AMR_ExAC TINYTEXT NOT NULL,
    Clinical_Significance_ClinVar TINYTEXT NOT NULL,
    AA_Allele_Frequency_ESP TINYTEXT NOT NULL,
    Protein_Change TINYTEXT NOT NULL,
    Variant_in_exLOVD TINYTEXT NOT NULL,
    EA_Allele_Frequency_ESP TINYTEXT NOT NULL,
    HGVS_RNA TINYTEXT NOT NULL,
    Clinical_significance_citations_ENIGMA TINYTEXT NOT NULL,
    Variant_effect_LOVD TINYTEXT NOT NULL,
    Polyphen_Prediction TINYTEXT NOT NULL,
    Data_Release_id TINYTEXT NOT NULL,
    Hg37_Start TINYTEXT NOT NULL,
    Hg36_Start TINYTEXT NOT NULL,
    Sift_Score TINYTEXT NOT NULL,
    Genomic_Coordinate_hg38 TINYTEXT NOT NULL,
    Alt TINYTEXT NOT NULL,
    Literature_citation_BIC TEXT NOT NULL,
    Variant_haplotype_LOVD TINYTEXT NOT NULL,
    Allele_frequency_NFE_ExAC TINYTEXT NOT NULL,
    Hg38_Start TINYTEXT NOT NULL,
    Pos TINYTEXT NOT NULL,
    Date_last_evaluated_ENIGMA TINYTEXT NOT NULL,
    Allele_number_SAS_ExAC TINYTEXT NOT NULL,
    Allele_number_AMR_ExAC TINYTEXT NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);
