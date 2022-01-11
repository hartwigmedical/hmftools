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
DROP TABLE IF EXISTS oncokbConsequenceBiological;
DROP TABLE IF EXISTS oncokbGeneBiological;
DROP TABLE IF EXISTS oncokbGeneAliasBiological;
DROP TABLE IF EXISTS oncokbClinical;
DROP TABLE IF EXISTS oncokbDrugAbstractClinical;
DROP TABLE IF EXISTS oncokbVariantClinical;
DROP TABLE IF EXISTS oncokbConsequenceClinical;
DROP TABLE IF EXISTS oncokbGeneClinical;
DROP TABLE IF EXISTS oncokbGeneAliasClinical;
DROP TABLE IF EXISTS civic;
DROP TABLE IF EXISTS civicAssertion;
DROP TABLE IF EXISTS civicHGVSExpression;
DROP TABLE IF EXISTS civicClinVarEntry;
DROP TABLE IF EXISTS civicVariantAlias;
DROP TABLE IF EXISTS civicVariantType;
DROP TABLE IF EXISTS civicProvisionalValue;
DROP TABLE IF EXISTS civicCoordinates;
DROP TABLE IF EXISTS civicVariantGroup;
DROP TABLE IF EXISTS civicVariantGroupVariant;
DROP TABLE IF EXISTS civicVariantGroupCoordinates;
DROP TABLE IF EXISTS civicVariantGroupType;
DROP TABLE IF EXISTS civicEvidenceItem;
DROP TABLE IF EXISTS civicDrug;
DROP TABLE IF EXISTS civicDisease;
DROP TABLE IF EXISTS civicEvidenceItemSource;
DROP TABLE IF EXISTS civicEvidenceItemPublication;
DROP TABLE IF EXISTS civicEvidenceItemClinicalTrial;
DROP TABLE IF EXISTS civicSource;
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
DROP TABLE IF EXISTS molecularMatch;
DROP TABLE IF EXISTS molecularMatchAst;
DROP TABLE IF EXISTS molecularMatchAstLeft;
DROP TABLE IF EXISTS molecularMatchAstLeftLeft;
DROP TABLE IF EXISTS molecularMatchAstLeftRight;
DROP TABLE IF EXISTS molecularMatchAstRight;
DROP TABLE IF EXISTS molecularMatchAstRightLeft;
DROP TABLE IF EXISTS molecularMatchAstRightRight;
DROP TABLE IF EXISTS molecularMatchInstitution;
DROP TABLE IF EXISTS molecularMatchExternalId;
DROP TABLE IF EXISTS molecularMatchIncludeGene1;
DROP TABLE IF EXISTS molecularMatchIncludeFinding1;
DROP TABLE IF EXISTS molecularMatchIncludeCondition1;
DROP TABLE IF EXISTS molecularMatchIncludeMutation1;
DROP TABLE IF EXISTS molecularMatchIncludeDrug1;
DROP TABLE IF EXISTS molecularMatchIncludeDrugClass1;
DROP TABLE IF EXISTS molecularMatchIncludeResistance1;
DROP TABLE IF EXISTS molecularMatchIncludeStage0;
DROP TABLE IF EXISTS molecularMatchIncludeGene0;
DROP TABLE IF EXISTS molecularMatchIncludeCondition0;
DROP TABLE IF EXISTS molecularMatchIncludeMutation0;
DROP TABLE IF EXISTS molecularMatchCriteriaMet;
DROP TABLE IF EXISTS molecularMatchSource;
DROP TABLE IF EXISTS molecularMatchTierExplanation;
DROP TABLE IF EXISTS molecularMatchTherapeuticContext;
DROP TABLE IF EXISTS molecularMatchTag;
DROP TABLE IF EXISTS molecularMatchVariantInfo;
DROP TABLE IF EXISTS molecularMatchVariantInfoConsequence;
DROP TABLE IF EXISTS molecularMatchVariantInfoFusion;
DROP TABLE IF EXISTS molecularMatchVariantInfoLocation;
DROP TABLE IF EXISTS molecularMatchVariantInfoLocationExonNumber;
DROP TABLE IF EXISTS molecularMatchClassification;
DROP TABLE IF EXISTS molecularMatchClassificationTranscript;
DROP TABLE IF EXISTS molecularMatchClassificationChromosome;
DROP TABLE IF EXISTS molecularMatchClassificationStart;
DROP TABLE IF EXISTS molecularMatchClassificationEnd;
DROP TABLE IF EXISTS molecularMatchClassificationRef;
DROP TABLE IF EXISTS molecularMatchClassificationAlt;
DROP TABLE IF EXISTS molecularMatchClassificationNucleotideChange;
DROP TABLE IF EXISTS molecularMatchClassificationExon;
DROP TABLE IF EXISTS molecularMatchClassificationExonicFunc;
DROP TABLE IF EXISTS molecularMatchClassificationPathology;
DROP TABLE IF EXISTS molecularMatchClassificationSource;
DROP TABLE IF EXISTS molecularMatchClassificationCosmicId;
DROP TABLE IF EXISTS molecularMatchClassificationDbSNP;
DROP TABLE IF EXISTS molecularMatchClassificationPopFreqMax;
DROP TABLE IF EXISTS molecularMatchClassificationParent;
DROP TABLE IF EXISTS molecularMatchClassificationParentTranscript;
DROP TABLE IF EXISTS molecularMatchCriteriaUnmet;
DROP TABLE IF EXISTS molecularMatchPrevalence;
DROP TABLE IF EXISTS molecularMatchMutation;
DROP TABLE IF EXISTS molecularMatchMutationMutationType;
DROP TABLE IF EXISTS molecularMatchMutationSource;
DROP TABLE IF EXISTS molecularMatchMutationSynonym;
DROP TABLE IF EXISTS molecularMatchMutationPathology;
DROP TABLE IF EXISTS molecularMatchMutationCDNA;
DROP TABLE IF EXISTS molecularMatchMutationTranscriptConsequence;
DROP TABLE IF EXISTS molecularMatchMutationTranscriptConsequenceExonNumber;
DROP TABLE IF EXISTS molecularMatchMutationParent;
DROP TABLE IF EXISTS molecularMatchMutationParentTranscript;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocation;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationGene;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationFullAA;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationClinVarDisease;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationClinVarSig;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationClinVarStatus;
DROP TABLE IF EXISTS molecularMatchMutationWGSALocationClinVarDbId;
DROP TABLE IF EXISTS molecularMatchMutationWGSAMap;
DROP TABLE IF EXISTS molecularMatchMutationWGSAMapSynonym;
DROP TABLE IF EXISTS molecularMatchMutationWGSAMapProtCoord;
DROP TABLE IF EXISTS molecularMatchMutationGRCh37Loc;
DROP TABLE IF EXISTS molecularMatchMutationGRCh37LocConsequence;
DROP TABLE IF EXISTS molecularMatchMutationGRCh37LocConsequenceTxSite;
DROP TABLE IF EXISTS molecularMatchMutationGRCh37LocConsequenceExonNumber;
DROP TABLE IF EXISTS molecularMatchMutationFusion;
DROP TABLE IF EXISTS molecularMatchMutationFusionAChromosome;
DROP TABLE IF EXISTS molecularMatchMutationFusionABand;
DROP TABLE IF EXISTS molecularMatchMutationFusionAGene;
DROP TABLE IF EXISTS molecularMatchMutationFusionACoord;
DROP TABLE IF EXISTS molecularMatchMutationFusionATranscript;
DROP TABLE IF EXISTS molecularMatchMutationFusionAOrientation;
DROP TABLE IF EXISTS molecularMatchMutationFusionAGenomicRegion;
DROP TABLE IF EXISTS molecularMatchMutationFusionBChromosome;
DROP TABLE IF EXISTS molecularMatchMutationFusionBBand;
DROP TABLE IF EXISTS molecularMatchMutationFusionBGene;
DROP TABLE IF EXISTS molecularMatchMutationFusionBCoord;
DROP TABLE IF EXISTS molecularMatchMutationFusionBTranscript;
DROP TABLE IF EXISTS molecularMatchMutationFusionBOrientation;
DROP TABLE IF EXISTS molecularMatchMutationFusionBGenomicRegion;
DROP TABLE IF EXISTS molecularMatchMutationFusionInsert;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfo;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon1;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon2;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon3;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon4;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon5;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon6;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon7;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon8;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon9;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon10;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon11;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon12;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon13;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon14;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon15;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon16;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon17;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon18;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon19;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon20;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon21;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon22;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon23;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon24;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon25;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon26;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon27;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon28;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon29;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon30;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon31;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon32;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon33;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon34;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon35;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon36;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon37;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon38;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon39;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon40;
DROP TABLE IF EXISTS molecularMatchMutationExonsInfoBoundaryExon41;
DROP TABLE IF EXISTS pmkb;
DROP TABLE IF EXISTS pmkbTissue;
DROP TABLE IF EXISTS pmkbTumor;
DROP TABLE IF EXISTS pmkbVariant;
DROP TABLE IF EXISTS pmkbGene;
DROP TABLE IF EXISTS molecularMatchTrials;
DROP TABLE IF EXISTS molecularMatchTrialsAlteration;
DROP TABLE IF EXISTS molecularMatchTrialsIntervention;
DROP TABLE IF EXISTS molecularMatchTrialsOtherName;
DROP TABLE IF EXISTS molecularMatchTrialsArmGroupLabel;
DROP TABLE IF EXISTS molecularMatchTrialsOverallContact;
DROP TABLE IF EXISTS molecularMatchTrialsTag;
DROP TABLE IF EXISTS molecularMatchTrialsLocation;
DROP TABLE IF EXISTS molecularMatchTrialsContact;
DROP TABLE IF EXISTS molecularMatchTrialsGeo;
DROP TABLE IF EXISTS molecularMatchTrialsSubLocation;
DROP TABLE IF EXISTS molecularMatchTrialsCoordinates;
DROP TABLE IF EXISTS jaxTrials;
DROP TABLE IF EXISTS jaxTrialsMolecularProfile;
DROP TABLE IF EXISTS jaxTrialsIndication;
DROP TABLE IF EXISTS jaxTrialsTherapy;
DROP TABLE IF EXISTS jax;
DROP TABLE IF EXISTS jaxMolecularProfile;
DROP TABLE IF EXISTS jaxIndication;
DROP TABLE IF EXISTS jaxTherapy;
DROP TABLE IF EXISTS jaxReference;
DROP TABLE IF EXISTS cgi;
DROP TABLE IF EXISTS cgiTranscript;
DROP TABLE IF EXISTS cgiIndividualMutation;
DROP TABLE IF EXISTS cgiGDNA;
DROP TABLE IF EXISTS cgiCDNA;
DROP TABLE IF EXISTS cgiInfo;
DROP TABLE IF EXISTS cgiRegion;
DROP TABLE IF EXISTS cgiStrand;
DROP TABLE IF EXISTS brca;
DROP TABLE IF EXISTS brcaAnnotation1000Genomes;
DROP TABLE IF EXISTS brcaAnnotationBIC;
DROP TABLE IF EXISTS brcaAnnotationClinVar;
DROP TABLE IF EXISTS brcaAnnotationENIGMA;
DROP TABLE IF EXISTS brcaAnnotationESP;
DROP TABLE IF EXISTS brcaAnnotationExAC;
DROP TABLE IF EXISTS brcaAnnotationExLOVD;
DROP TABLE IF EXISTS brcaAnnotationLOVD;

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
    gene varchar(255) NOT NULL,
    entrezId varchar(255) NOT NULL,
    clinicalManifestation varchar(255) NOT NULL,
    responseType varchar(255) NOT NULL,
    evidenceLabel varchar(255) NOT NULL,
    drugLabels varchar(255) NOT NULL,
    germlineOrSomatic varchar(255) NOT NULL,
    publicationUrl varchar(255) NOT NULL,
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
    oncokbId int NOT NULL,
    gene varchar(255) NOT NULL,
    entrezGeneId varchar(255) NOT NULL,
    isoform varchar(255) NOT NULL,
    refSeq varchar(255) NOT NULL,
    oncogenic varchar(255) NOT NULL,
    mutationEffect varchar(255) NOT NULL,
    mutationEffectPmids varchar(255) NOT NULL,
    mutationEffectAbstracts varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbId) REFERENCES oncokb(id)
);

CREATE TABLE oncokbVariantBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbBiologicalId int NOT NULL,
    name varchar(255) NOT NULL,
    alteration varchar(255) NOT NULL,
    proteinStart varchar(255) NOT NULL,
    proteinEnd varchar(255) NOT NULL,
    refResidues varchar(255),
    variantResidues varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbBiologicalId) REFERENCES oncokbBiological(id)
);

CREATE TABLE oncokbConsequenceBiological
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
    oncokbVariantBiologicalId int NOT NULL,
    hugoSymbol varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    entrezGeneId varchar(255) NOT NULL,
    curatedIsoform varchar(255),
    curatedRefSeq varchar(255),
    oncogene varchar(255) NOT NULL,
    tsg varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbVariantBiologicalId) REFERENCES oncokbVariantBiological(id)
);

CREATE TABLE oncokbGeneAliasBiological
(   id int NOT NULL AUTO_INCREMENT,
    oncokbGeneBiologicalId int NOT NULL,
    geneAlias varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbGeneBiologicalId) REFERENCES oncokbGeneBiological(id)
);

CREATE TABLE oncokbClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbId int NOT NULL,
    gene varchar(255) NOT NULL,
    entrezGeneId varchar(255) NOT NULL,
    isoform varchar(255) NOT NULL,
    refSeq varchar(255) NOT NULL,
    cancerType varchar(255) NOT NULL,
    drug varchar(255) NOT NULL,
    drugPmids varchar(255) NOT NULL,
    level varchar(255) NOT NULL,
    levelLabel varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbId) REFERENCES oncokb(id)
);

CREATE TABLE oncokbDrugAbstractClinical
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
    name varchar(255) NOT NULL,
    alteration varchar(255) NOT NULL,
    proteinStart varchar(255) NOT NULL,
    proteinEnd varchar(255) NOT NULL,
    refResidues varchar(255),
    variantResidues varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbClinicalId) REFERENCES oncokbClinical(id)
);

CREATE TABLE oncokbConsequenceClinical
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
    oncokbVariantClinicalId int NOT NULL,
    hugoSymbol varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    entrezGeneId varchar(255) NOT NULL,
    curatedIsoform varchar(255),
    curatedRefSeq varchar(255),
    oncogene varchar(255) NOT NULL,
    tsg varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbVariantClinicalId) REFERENCES oncokbVariantClinical(id)
);

CREATE TABLE oncokbGeneAliasClinical
(   id int NOT NULL AUTO_INCREMENT,
    oncokbGeneClinicalId int NOT NULL,
    geneAlias varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (oncokbGeneClinicalId) REFERENCES oncokbGeneClinical(id)
);

CREATE TABLE civic
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    entrezId varchar(255) NOT NULL,
    entrezName varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    civicActionabilityScore varchar(255),
    alleleRegistryId varchar(255),
    idCivic varchar(255) NOT NULL,
    geneId varchar(255) NOT NULL,
    description varchar(2000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE civicAssertion
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    assertion varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicHGVSExpression
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    hgvsExpression varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicClinVarEntry
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    clinVarEntry varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantAlias
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    variantAlias varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantType
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    description varchar(255) NOT NULL,
    url varchar(255) NOT NULL,
    soId varchar(255) NOT NULL,
    idVariantType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicProvisionalValue
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
    chromosome varchar(255),
    start varchar(255),
    stop varchar(255),
    referenceBases varchar(255),
    variantBases varchar(255),
    representativeTranscript varchar(255),
    ensemblVersion varchar(255),
    referenceBuild varchar(255),
    chromosome2 varchar(255),
    start2 varchar(255),
    stop2 varchar(255),
    representativeTranscript2 varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantGroup
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    idVariantGroup varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicVariantGroupVariant
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantGroupId int NOT NULL,
    entrezId varchar(255) NOT NULL,
    entrezName varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    civicActionabilityScore varchar(255),
    idVariant varchar(255) NOT NULL,
    geneId varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantGroupId) REFERENCES civicVariantGroup(id)
);

CREATE TABLE civicVariantGroupCoordinates
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantGroupVariantId int NOT NULL,
    chromosome varchar(255),
    start varchar(255),
    stop varchar(255),
    referenceBases varchar(255),
    variantBases varchar(255),
    representativeTranscript varchar(255),
    ensemblVersion varchar(255),
    referenceBuild varchar(255),
    chromosome2 varchar(255),
    start2 varchar(255),
    stop2 varchar(255),
    representativeTranscript2 varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantGroupVariantId) REFERENCES civicVariantGroupVariant(id)
);

CREATE TABLE civicVariantGroupType
(   id int NOT NULL AUTO_INCREMENT,
    civicVariantGroupVariantId int NOT NULL,
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    url varchar(255) NOT NULL,
    soId varchar(255) NOT NULL,
    idVariantType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicVariantGroupVariantId) REFERENCES civicVariantGroupVariant(id)
);

CREATE TABLE civicEvidenceItem
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    status varchar(255) NOT NULL,
    rating varchar(255),
    evidenceType varchar(255) NOT NULL,
    evidenceLevel varchar(255) NOT NULL,
    evidenceDirection varchar(255),
    drugInteractionType varchar(255),
    variantOrigin varchar(255),
    clinicalSignificance varchar(255),
    openChangeCount varchar(255) NOT NULL,
    description varchar(1500) NOT NULL,
    variantId varchar(255),
    idEvidenceItem varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicDrug
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemId int NOT NULL,
    name varchar(255) NOT NULL,
    pubchemId varchar(255),
    idDrug varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemId) REFERENCES civicEvidenceItem(id)
);

CREATE TABLE civicDisease
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemId int NOT NULL,
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    doid varchar(255),
    url varchar(255),
    idDisease varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemId) REFERENCES civicEvidenceItem(id)
);

CREATE TABLE civicEvidenceItemSource
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemId int NOT NULL,
    name varchar(1000),
    status varchar(255) NOT NULL,
    openAccess varchar(255),
    journal varchar(255),
    fullJournalTitle varchar(255),
    citation varchar(255) NOT NULL,
    pmcId varchar(255),
    sourceUrl varchar(255) NOT NULL,
    pubmedId varchar(255) NOT NULL,
    isReview varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemId) REFERENCES civicEvidenceItem(id)
);

CREATE TABLE civicEvidenceItemPublication
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemSourceId int NOT NULL,
    year varchar(255),
    month varchar(255),
    day varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemSourceId) REFERENCES civicEvidenceItemSource(id)
);

CREATE TABLE civicEvidenceItemClinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    civicEvidenceItemSourceId int NOT NULL,
    name varchar(255),
    nctId varchar(255),
    clinicalTrialUrl varchar(255),
    description varchar(2000),
    PRIMARY KEY (id),
    FOREIGN KEY (civicEvidenceItemSourceId) REFERENCES civicEvidenceItemSource(id)
);

CREATE TABLE civicSource
(   id int NOT NULL AUTO_INCREMENT,
    civicId int NOT NULL,
    name varchar(255),
    status varchar(255) NOT NULL,
    openAccess varchar(255),
    journal varchar(255),
    fullJournalTitle varchar(255),
    citation varchar(255) NOT NULL,
    pmcId varchar(255),
    sourceUrl varchar(255) NOT NULL,
    pubmedId varchar(255) NOT NULL,
    isReview varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicId) REFERENCES civic(id)
);

CREATE TABLE civicPublication
(   id int NOT NULL AUTO_INCREMENT,
    civicSourceId int NOT NULL,
    year varchar(255),
    month varchar(255),
    day varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (civicSourceId) REFERENCES civicSource(id)
);

CREATE TABLE civicClinicalTrial
(   id int NOT NULL AUTO_INCREMENT,
    civicSourceId int NOT NULL,
    name varchar(255),
    nctId varchar(255),
    clinicalTrialUrl varchar(255),
    description varchar(1000),
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
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    affiliation varchar(255),
    featuredExpert varchar(255) NOT NULL,
    areaOfExpertise varchar(255),
    bio varchar(1000),
    url varchar(255),
    createdAt varchar(255) NOT NULL,
    lastSeenAt varchar(255),
    avatarUrl varchar(255) NOT NULL,
    twitterHandle varchar(255),
    facebookProfile varchar(255),
    linkedinProfile varchar(255),
    orcid varchar(255),
    signupComplete varchar(255),
    acceptedLicense varchar(255),
    idUser varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnId) REFERENCES civicLastCommentedOn(id)
);

CREATE TABLE civicLastCommentedOnAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnUserId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnUserId) REFERENCES civicLastCommentedOnUser(id)
);

CREATE TABLE civicLastCommentedOnOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnUserId int NOT NULL,
    name varchar(255) NOT NULL,
    url varchar(255) NOT NULL,
    idOrganization varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastCommentedOnUserId) REFERENCES civicLastCommentedOnUser(id)
);

CREATE TABLE civicLastCommentedOnProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastCommentedOnOrganizationId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
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
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    affiliation varchar(255),
    featuredExpert varchar(255) NOT NULL,
    areaOfExpertise varchar(255),
    bio varchar(1500),
    url varchar(255),
    createdAt varchar(255) NOT NULL,
    lastSeenAt varchar(255),
    avatarUrl varchar(255) NOT NULL,
    twitterHandle varchar(255),
    facebookProfile varchar(255),
    linkedinProfile varchar(255),
    orcid varchar(255),
    signupComplete varchar(255),
    acceptedLicense varchar(255),
    idUser varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedId) REFERENCES civicLastModified(id)
);

CREATE TABLE civicLastModifiedAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedId) REFERENCES civicLastModified(id)
);

CREATE TABLE civicLastModifiedOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedUserId int NOT NULL,
    name varchar(255),
    url varchar(255),
    idOrganization varchar(255),
    description varchar(1000),
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastModifiedUserId) REFERENCES civicLastModifiedUser(id)
);

CREATE TABLE civicLastModifiedProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastModifiedOrganizationId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
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
    name varchar(255) NOT NULL,
    displayName varchar(255) NOT NULL,
    role varchar(255) NOT NULL,
    affiliation varchar(255),
    featuredExpert varchar(255) NOT NULL,
    areaOfExpertise varchar(255),
    bio varchar(1000),
    url varchar(255),
    createdAt varchar(255) NOT NULL,
    lastSeenAt varchar(255),
    avatarUrl varchar(255) NOT NULL,
    twitterHandle varchar(255),
    facebookProfile varchar(255),
    linkedinProfile varchar(255),
    orcid varchar(255),
    signupComplete varchar(255),
    acceptedLicense varchar(255),
    idUser varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedId) REFERENCES civicLastReviewed(id)
);

CREATE TABLE civicLastReviewedAvatars
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedId) REFERENCES civicLastReviewed(id)
);

CREATE TABLE civicLastReviewedOrganization
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedUserId int NOT NULL,
    name varchar(255) NOT NULL,
    url varchar(255) NOT NULL,
    idOrganization varchar(255) NOT NULL,
    description varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedUserId) REFERENCES civicLastReviewedUser(id)
);

CREATE TABLE civicLastReviewedProfileImage
(   id int NOT NULL AUTO_INCREMENT,
    civicLastReviewedOrganizationId int NOT NULL,
    x14 varchar(255) NOT NULL,
    x32 varchar(255) NOT NULL,
    x64 varchar(255) NOT NULL,
    x128 varchar(255) NOT NULL,
    x256 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (civicLastReviewedOrganizationId) REFERENCES civicLastReviewedOrganization(id)
);

CREATE TABLE molecularMatch
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    direction varchar(255) NOT NULL,
    biomarkerClass varchar(255) NOT NULL,
    score varchar(255) NOT NULL,
    clinicalSignificance varchar(255) NOT NULL,
    tier varchar(255) NOT NULL,
    ampcap varchar(255) NOT NULL,
    civicValue varchar(255) NOT NULL,
    regulatoryBody varchar(255) NOT NULL,
    regulatoryBodyApproved varchar(255) NOT NULL,
    guidelineBody varchar(255),
    guidelineVersion varchar(255),
    noTherapyAvailable varchar(255),
    sixtier varchar(255) NOT NULL,
    mvld varchar(255) NOT NULL,
    autoGenerateNarrative varchar(255) NOT NULL,
    narrative varchar(1000) NOT NULL,
    expression varchar(1000) NOT NULL,
    customer varchar(255) NOT NULL,
    version varchar(255) NOT NULL,
    idMolecularMatch varchar(255) NOT NULL,
    uniqueKey varchar(255) NOT NULL,
    hashKey varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE molecularMatchAst
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchAstLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstId) REFERENCES molecularMatchAst(id)
);

CREATE TABLE molecularMatchAstLeftLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstLeftId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstLeftId) REFERENCES molecularMatchAstLeft(id)
);

CREATE TABLE molecularMatchAstLeftRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstLeftId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstLeftId) REFERENCES molecularMatchAstLeft(id)
);

CREATE TABLE molecularMatchAstRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(1000),
    value varchar(1000),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstId) REFERENCES molecularMatchAst(id)
);

CREATE TABLE molecularMatchAstRightLeft
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstRightId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstRightId) REFERENCES molecularMatchAstRight(id)
);

CREATE TABLE molecularMatchAstRightRight
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchAstRightId int NOT NULL,
    type varchar(255) NOT NULL,
    raw varchar(255),
    value varchar(255),
    operator varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchAstRightId) REFERENCES molecularMatchAstRight(id)
);

CREATE TABLE molecularMatchInstitution
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    institution varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchExternalId
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    externalId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeGene1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeGene1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeFinding1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeFinding1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeCondition1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeCondition1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeMutation1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeMutation1 varchar(1000) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeDrug1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeDrug1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeDrugClass1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeDrugClass1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeResistance1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeResistance1 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeStage0
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeStage0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeGene0
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeGene0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeCondition0
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeCondition0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchIncludeMutation0
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    includeMutation0 varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchCriteriaMet
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    criteriaMet varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchSource
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    subType varchar(255),
    valid varchar(255) NOT NULL,
    pubId varchar(255) NOT NULL,
    link varchar(255) NOT NULL,
    trialId varchar(255),
    trialPhase varchar(255),
    year varchar(255) NOT NULL,
    functionalConsequence varchar(255),
    institution varchar(255),
    trustRating varchar(255),
    suppress varchar(255) NOT NULL,
    idSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchTierExplanation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    tier varchar(255) NOT NULL,
    step varchar(255) NOT NULL,
    message varchar(1500) NOT NULL,
    success varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchTherapeuticContext
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    name varchar(255) NOT NULL,
    facet varchar(255) NOT NULL,
    suppress varchar(255),
    valid varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchTag
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    term varchar(1000) NOT NULL,
    facet varchar(255) NOT NULL,
    filterType varchar(255),
    priority varchar(255) NOT NULL,
    transcript varchar(255),
    valid varchar(255),
    generatedBy varchar(1000),
    generatedByTerm varchar(1000),
    isNew varchar(255),
    primaryValue varchar(255),
    custom varchar(255),
    suppress varchar(255),
    manualSuppress varchar(255),
    composite varchar(255),
    compositeKey varchar(3000),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchVariantInfo
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    name varchar(1000) NOT NULL,
    gene varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    classification varchar(255) NOT NULL,
    geneFusionPartner varchar(255) NOT NULL,
    cosmicId varchar(255),
    popFreqMax varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchVariantInfoConsequence
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    consequence varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoFusion
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    chr varchar(255) NOT NULL,
    referenceGenome varchar(255) NOT NULL,
    LBPWREP varchar(255) NOT NULL,
    LBPWLEP varchar(255) NOT NULL,
    RBPWREP varchar(255) NOT NULL,
    RBPWLEP varchar(255) NOT NULL,
    intronNumber varchar(255) NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoLocation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoId int NOT NULL,
    chr varchar(255) NOT NULL,
    start varchar(255) NOT NULL,
    stop varchar(255) NOT NULL,
    ref varchar(1000),
    alt varchar(255),
    cdna varchar(255),
    aminoAcidChange varchar(1000),
    referenceGenome varchar(255),
    strand varchar(255),
    intronNumber varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoId) REFERENCES molecularMatchVariantInfo(id)
);

CREATE TABLE molecularMatchVariantInfoLocationExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchVariantInfoLocationId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchVariantInfoLocationId) REFERENCES molecularMatchVariantInfoLocation(id)
);

CREATE TABLE molecularMatchClassification
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    name varchar(1000) NOT NULL,
    geneSymbol varchar(255),
    expandGeneSearch varchar(255),
    transcript varchar(255),
    classification varchar(255) NOT NULL,
    classificationOverride varchar(255),
    copyNumberType varchar(255),
    drugsApprovedOnLabelCount varchar(255),
    drugsApprovedOffLabelCount varchar(255),
    drugsExperimentalCount varchar(255),
    trialCount varchar(255),
    publicationCount varchar(255),
    rootTerm varchar(1000),
    alias varchar(1000),
    priority varchar(255),
    description varchar(1500),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchClassificationTranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationChromosome
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    chromosome varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationStart
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    start varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationEnd
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    end varchar(255) NOT NULL,
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

CREATE TABLE molecularMatchClassificationAlt
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    alt varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationNucleotideChange
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    nucleotideChange varchar(255) NOT NULL,
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

CREATE TABLE molecularMatchClassificationExonicFunc
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    exonicFunc varchar(255) NOT NULL,
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

CREATE TABLE molecularMatchClassificationSource
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationCosmicId
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    cosmicId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationDbSNP
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

CREATE TABLE molecularMatchClassificationParent
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationId int NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255),
    actionableParent varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationId) REFERENCES molecularMatchClassification(id)
);

CREATE TABLE molecularMatchClassificationParentTranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchClassificationParentId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchClassificationParentId) REFERENCES molecularMatchClassificationParent(id)
);

CREATE TABLE molecularMatchCriteriaUnmet
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    term varchar(1000) NOT NULL,
    filterType varchar(255) NOT NULL,
    priority varchar(255) NOT NULL,
    facet varchar(255) NOT NULL,
    valid varchar(255),
    transcript varchar(255),
    isNew varchar(255),
    generatedBy varchar(255),
    generatedByTerm varchar(255),
    suppress varchar(255) NOT NULL,
    manualSuppress varchar(255),
    primaryValue varchar(255),
    compositeKey varchar(3000) NOT NULL,
    custom varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchPrevalence
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    studyId varchar(255) NOT NULL,
    count varchar(255) NOT NULL,
    samples varchar(255) NOT NULL,
    percent varchar(255) NOT NULL,
    molecular varchar(255),
    conditionValue varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchMutation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchId int NOT NULL,
    geneSymbol varchar(255) NOT NULL,
    name varchar(1000) NOT NULL,
    transcriptRecognized varchar(255),
    transcript varchar(255),
    longestTranscript varchar(255),
    uniprotTranscript varchar(255),
    description varchar(1500) NOT NULL,
    src varchar(255) NOT NULL,
    idMutation varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchId) REFERENCES molecularMatch(id)
);

CREATE TABLE molecularMatchMutationTranscriptConsequence
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    chr varchar(255),
    start varchar(255),
    stop varchar(255),
    ref varchar(1000),
    alt varchar(255),
    referenceGenome varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    strand varchar(255) NOT NULL,
    cdna varchar(255),
    aminoAcidChange varchar(255),
    intronNumber varchar(255),
    suppress varchar(255) NOT NULL,
    custom varchar(255) NOT NULL,
    validated varchar(255) NOT NULL,
    compositeKey varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationTranscriptConsequenceExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationTranscriptConsequenceId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationTranscriptConsequenceId) REFERENCES molecularMatchMutationTranscriptConsequence(id)
);

CREATE TABLE molecularMatchMutationParent
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    name varchar(255) NOT NULL,
    type varchar(255),
    actionableParent varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationParentTranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationParentId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationParentId) REFERENCES molecularMatchMutationParent(id)
);

CREATE TABLE molecularMatchMutationMutationType
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    mutationType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationSource
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    source varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationSynonym
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    synonym varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationPathology
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    pathology varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationCDNA
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    cDNA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationWGSALocation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    chr varchar(255) NOT NULL,
    start varchar(255) NOT NULL,
    end varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    chrStartRefAlt varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    nucleotideChange varchar(255) NOT NULL,
    aa varchar(255),
    exonicFunc varchar(255),
    popFreqMax varchar(255) NOT NULL,
    exacAFR varchar(255),
    exacAMR varchar(255),
    exacEAS varchar(255),
    exacFIN varchar(255),
    exacNFE varchar(255),
    exacSAS varchar(255),
    exacFreq varchar(255),
    g1000AFR varchar(255),
    g1000AMR varchar(255),
    g1000EAS varchar(255),
    g1000EUR varchar(255),
    g1000SAS varchar(255),
    g1000ALL varchar(255),
    fathmm varchar(255) NOT NULL,
    fathmmPred varchar(255) NOT NULL,
    esp6500siAA varchar(255),
    esp6500siEA varchar(255),
    dbSNP varchar(255),
    cosmicId varchar(255),
    phyloP46wayPlacental varchar(255) NOT NULL,
    phyloP100wayVertebrate varchar(255) NOT NULL,
    siPhy29wayLogOdds varchar(255) NOT NULL,
    gwasSNP varchar(255),
    gwasDIS varchar(255),
    gwasPubmed varchar(255),
    gerpRS varchar(255) NOT NULL,
    func varchar(255) NOT NULL,
    wgRna varchar(255),
    targetScanS varchar(255),
    keyValue varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationGene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationFullAA
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    fullAA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationClinVarDisease
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    clinVarDisease varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationClinVarSig
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    clinVarSig varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationClinVarStatus
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    clinVarStatus varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSALocationClinVarDbId
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSALocationId int NOT NULL,
    clinVarDbId varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSALocationId) REFERENCES molecularMatchMutationWGSALocation(id)
);

CREATE TABLE molecularMatchMutationWGSAMap
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    name varchar(255) NOT NULL,
    gene varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    exon varchar(255),
    grch37ChrStartRefAlt varchar(255) NOT NULL,
    nucleotideChange varchar(255) NOT NULL,
    aa varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationWGSAMapSynonym
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSAMapId int NOT NULL,
    synonym varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSAMapId) REFERENCES molecularMatchMutationWGSAMap(id)
);

CREATE TABLE molecularMatchMutationWGSAMapProtCoord
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationWGSAMapId int NOT NULL,
    protCoord varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationWGSAMapId) REFERENCES molecularMatchMutationWGSAMap(id)
);

CREATE TABLE molecularMatchMutationGRCh37Loc
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    chr varchar(255),
    start varchar(255),
    stop varchar(255),
    ref varchar(1000),
    alt varchar(255),
    strand varchar(255) NOT NULL,
    validated varchar(255) NOT NULL,
    compositeKey varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationGRCh37LocConsequence
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationGRCh37LocId int NOT NULL,
    transcript varchar(255) NOT NULL,
    cdna varchar(255),
    aminoAcidChange varchar(255),
    intronNumber varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationGRCh37LocId) REFERENCES molecularMatchMutationGRCh37Loc(id)
);

CREATE TABLE molecularMatchMutationGRCh37LocConsequenceTxSite
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationGRCh37LocConsequenceId int NOT NULL,
    txSite varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationGRCh37LocConsequenceId) REFERENCES molecularMatchMutationGRCh37LocConsequence(id)
);

CREATE TABLE molecularMatchMutationGRCh37LocConsequenceExonNumber
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationGRCh37LocConsequenceId int NOT NULL,
    exonNumber varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationGRCh37LocConsequenceId) REFERENCES molecularMatchMutationGRCh37LocConsequence(id)
);

CREATE TABLE molecularMatchMutationFusion
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    source varchar(255),
    synonym varchar(255),
    paper varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationFusionAChromosome
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    chromosome varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionABand
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    band varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionAGene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionACoord
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    coord varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionATranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionAOrientation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    orientation varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionAGenomicRegion
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    num varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBChromosome
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    chromosome varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBBand
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    band varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBGene
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    gene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBCoord
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    coord varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBTranscript
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    transcript varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBOrientation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    orientation varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionBGenomicRegion
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    num varchar(255) NOT NULL,
    type varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationFusionInsert
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationFusionId int NOT NULL,
    ins varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationFusionId) REFERENCES molecularMatchMutationFusion(id)
);

CREATE TABLE molecularMatchMutationExonsInfo
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationId int NOT NULL,
    chr varchar(255) NOT NULL,
    transcript varchar(255) NOT NULL,
    txStart varchar(255),
    txEnd varchar(255),
    cdsStart varchar(255),
    cdsEnd varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationId) REFERENCES molecularMatchMutation(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon1
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon2
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon3
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon4
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon5
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon6
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon7
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon8
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon9
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon10
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon11
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon12
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon13
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon14
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon15
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon16
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon17
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon18
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon19
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon20
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon21
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon22
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon23
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon24
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon25
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon26
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon27
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon28
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon29
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon30
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon31
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon32
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon33
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon34
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon35
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon36
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon37
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon38
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon39
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon40
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE molecularMatchMutationExonsInfoBoundaryExon41
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchMutationExonsInfoId int NOT NULL,
    start varchar(255),
    end varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchMutationExonsInfoId) REFERENCES molecularMatchMutationExonsInfo(id)
);

CREATE TABLE pmkb
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE pmkbTissue
(   id int NOT NULL AUTO_INCREMENT,
    pmkbId int NOT NULL,
    name varchar(255) NOT NULL,
    idTissue varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (pmkbId) REFERENCES pmkb(id)
);

CREATE TABLE pmkbTumor
(   id int NOT NULL AUTO_INCREMENT,
    pmkbId int NOT NULL,
    name varchar(255) NOT NULL,
    idTumor varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (pmkbId) REFERENCES pmkb(id)
);

CREATE TABLE pmkbVariant
(   id int NOT NULL AUTO_INCREMENT,
    pmkbId int NOT NULL,
    name varchar(255),
    coordinates varchar(255),
    chromosome varchar(255),
    cytoband varchar(255),
    transcript varchar(255),
    effect varchar(255),
    codons varchar(255),
    exons varchar(255),
    dnaChange varchar(255),
    aminoAcidChange varchar(255),
    germline varchar(255),
    partnerGene varchar(255),
    cnvType varchar(255),
    chromosomeBasedCnv varchar(255),
    variantType varchar(255),
    cosmic varchar(255),
    description varchar(255),
    descriptionType varchar(255),
    notes varchar(255),
    idVariant varchar(255),
    PRIMARY KEY (id),
    FOREIGN KEY (pmkbId) REFERENCES pmkb(id)
);

CREATE TABLE pmkbGene
(   id int NOT NULL AUTO_INCREMENT,
    pmkbVariantId int NOT NULL,
    name varchar(255) NOT NULL,
    createdAt varchar(255) NOT NULL,
    updatedAt varchar(255) NOT NULL,
    activeInd varchar(255) NOT NULL,
    description varchar(255),
    externalId varchar(255) NOT NULL,
    idGene varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (pmkbVariantId) REFERENCES pmkbVariant(id)
);

CREATE TABLE molecularMatchTrials
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    status varchar(500) NOT NULL,
    startDate varchar(255),
    title varchar(1000) NOT NULL,
    briefTitle varchar(1000),
    studyType varchar(255) NOT NULL,
    score varchar(255) NOT NULL,
    link varchar(255) NOT NULL,
    phase varchar(255) NOT NULL,
    idMolecularMatchTrials varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE molecularMatchTrialsAlteration
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    molecularAlteration varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsIntervention
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    interventionName varchar(9000),
    interventionType varchar(250),
    description varchar(2000),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsOtherName
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsInterventionId int NOT NULL,
    otherName varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsInterventionId) REFERENCES molecularMatchTrialsIntervention(id)
);

CREATE TABLE molecularMatchTrialsArmGroupLabel
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsInterventionId int NOT NULL,
    armGroupLabel varchar(500),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsInterventionId) REFERENCES molecularMatchTrialsIntervention(id)
);

CREATE TABLE molecularMatchTrialsOverallContact
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    name varchar(250),
    type varchar(250),
    affiliation varchar(250),
    lastName varchar(250),
    email varchar(250),
    phone varchar(250),
    phoneExt varchar(250),
    street varchar(250),
    zip varchar(250),
    city varchar(250),
    country varchar(250),
    url varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsTag
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
    idTag varchar(250),
    manualPriority varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsLocation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsId int NOT NULL,
    status varchar(250) NOT NULL,
    name varchar(1000),
    lastName varchar(1000),
    email varchar(1000),
    phone varchar(250),
    phoneExt varchar(250),
    lastNameBackup varchar(250),
    emailBackup varchar(250),
    phoneBackup varchar(250),
    phoneExtBackup varchar(250),
    street varchar(250),
    city varchar(250),
    zip varchar(250),
    state varchar(250),
    country varchar(250),
    number varchar(250),
    poBox varchar(250),
    idLocation varchar(250),
    valid varchar(250),
    validMessage varchar(250),
    created varchar(250),
    lastUpdated varchar(250),
    failedGeocode varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsId) REFERENCES molecularMatchTrials(id)
);

CREATE TABLE molecularMatchTrialsContact
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationId int NOT NULL,
    name varchar(250),
    email varchar(250),
    phone varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationId) REFERENCES molecularMatchTrialsLocation(id)
);

CREATE TABLE molecularMatchTrialsGeo
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationId int NOT NULL,
    lat varchar(250),
    lon varchar(250),
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationId) REFERENCES molecularMatchTrialsLocation(id)
);

CREATE TABLE molecularMatchTrialsSubLocation
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsLocationId int NOT NULL,
    type varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsLocationId) REFERENCES molecularMatchTrialsLocation(id)
);

CREATE TABLE molecularMatchTrialsCoordinates
(   id int NOT NULL AUTO_INCREMENT,
    molecularMatchTrialsSubLocationId int NOT NULL,
    coordinates varchar(250) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (molecularMatchTrialsSubLocationId) REFERENCES molecularMatchTrialsSubLocation(id)
);

CREATE TABLE jaxTrials
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    nctId varchar(255) NOT NULL,
    title varchar(500) NOT NULL,
    variantRequirements varchar(255) NOT NULL,
    gender varchar(255),
    recruitment varchar(255) NOT NULL,
    phase varchar(255) NOT NULL,
    sponsors varchar(255) NOT NULL,
    updateDate varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE jaxTrialsMolecularProfile
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    requirementType varchar(255) NOT NULL,
    profileName varchar(255) NOT NULL,
    idMolecularProfile varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsId) REFERENCES jaxTrials(id)
);

CREATE TABLE jaxTrialsIndication
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    name varchar(255) NOT NULL,
    source varchar(255) NOT NULL,
    idIndication varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (jaxTrialsId) REFERENCES jaxTrials(id)
);

CREATE TABLE jaxTrialsTherapy
(   id int NOT NULL AUTO_INCREMENT,
    jaxTrialsId int NOT NULL,
    therapyName varchar(255) NOT NULL,
    idTherapy varchar(255) NOT NULL,
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

CREATE TABLE jaxIndication
(   id int NOT NULL AUTO_INCREMENT,
    jaxId int NOT NULL,
    source varchar(255) NOT NULL,
    idIndication varchar(255) NOT NULL,
    name varchar(255) NOT NULL,
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
    gene varchar(255) NOT NULL,
    biomarker varchar(1000) NOT NULL,
    alteration varchar(1000) NOT NULL,
    alterationType varchar(255) NOT NULL,
    association varchar(255) NOT NULL,
    drug varchar(255) NOT NULL,
    drugFamily varchar(255) NOT NULL,
    drugFullName varchar(255) NOT NULL,
    drugStatus varchar(255) NOT NULL,
    targeting varchar(255) NOT NULL,
    primaryTumorType varchar(255) NOT NULL,
    metastaticTumorType varchar(255) NOT NULL,
    evidenceLevel varchar(255) NOT NULL,
    source varchar(255) NOT NULL,
    curator varchar(255) NOT NULL,
    assayType varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE cgiTranscript
( id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    transcript varchar(255) NOT NULL,
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

CREATE TABLE cgiGDNA
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    gDNA varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE cgiCDNA
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    cDNA varchar(255) NOT NULL,
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

CREATE TABLE cgiStrand
(   id int NOT NULL AUTO_INCREMENT,
    cgiId int NOT NULL,
    strand varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (cgiId) REFERENCES cgi(id)
);

CREATE TABLE brca
(   id int NOT NULL AUTO_INCREMENT,
    viccEntryId int NOT NULL,
    geneSymbol varchar(255) NOT NULL,
    chr varchar(255) NOT NULL,
    pos varchar(255) NOT NULL,
    ref varchar(255) NOT NULL,
    alt varchar(255) NOT NULL,
    genomicCoordinateHg36 varchar(255) NOT NULL,
    hg36Start varchar(255) NOT NULL,
    hg36End varchar(255) NOT NULL,
    genomicCoordinateHg37 varchar(255) NOT NULL,
    hg37Start varchar(255) NOT NULL,
    hg37End varchar(255) NOT NULL,
    genomicCoordinateHg38 varchar(255) NOT NULL,
    hg38Start varchar(255) NOT NULL,
    hg38End varchar(255) NOT NULL,
    proteinChange varchar(255) NOT NULL,
    referenceSequence varchar(255) NOT NULL,
    synonyms varchar(1000) NOT NULL,
    hgvsCDNA varchar(255) NOT NULL,
    hgvsProtein varchar(255) NOT NULL,
    hgvsRNA varchar(255) NOT NULL,
    siftScore varchar(255) NOT NULL,
    siftPrediction varchar(255) NOT NULL,
    polyphenScore varchar(255) NOT NULL,
    polyphenPrediction varchar(255) NOT NULL,
    pathogenicityAll varchar(255) NOT NULL,
    pathogenicityExpert varchar(255) NOT NULL,
    alleleFrequency varchar(255) NOT NULL,
    maxAlleleFrequency varchar(255) NOT NULL,
    discordant varchar(255) NOT NULL,
    idBrca varchar(255) NOT NULL,
    changeTypeId varchar(255) NOT NULL,
    dataReleaseId varchar(255) NOT NULL,
    source varchar(255) NOT NULL,
    sourceURL varchar(2500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (viccEntryId) REFERENCES viccEntry(id)
);

CREATE TABLE brcaAnnotation1000Genomes (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantIn1000Genomes varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    alleleFrequency varchar(255) NOT NULL,
    afrAlleleFrequency varchar(255) NOT NULL,
    amrAlleleFrequency varchar(255) NOT NULL,
    easAlleleFrequency varchar(255) NOT NULL,
    eurAlleleFrequency varchar(255) NOT NULL,
    sasAlleleFrequency varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationBIC (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInBIC varchar(10) NOT NULL,
    bxId varchar(12500) NOT NULL,
    mutationType varchar(10) NOT NULL,
    clinicalClassification varchar(20) NOT NULL,
    clinicalImportance varchar(20) NOT NULL,
    nomenclature varchar(100) NOT NULL,
    ethnicity varchar(2000) NOT NULL,
    patientNationality varchar(400) NOT NULL,
    germlineOrSomatic varchar(20) NOT NULL,
    numberOfFamilyMemberCarryingMutation varchar(25) NOT NULL,
    literatureCitation varchar(600) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationClinVar (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInClinVar varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    clinicalSignificance varchar(255) NOT NULL,
    submitter varchar(2500) NOT NULL,
    method varchar(255) NOT NULL,
    alleleOrigin varchar(255) NOT NULL,
    scv varchar(500) NOT NULL,
    dateLastUpdated varchar(500) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationENIGMA (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInENIGMA varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    alleleOrigin varchar(255) NOT NULL,
    clinVarAccession varchar(255) NOT NULL,
    assertionMethod varchar(255) NOT NULL,
    assertionMethodCitation varchar(255) NOT NULL,
    collectionMethod varchar(255) NOT NULL,
    conditionCategory varchar(255) NOT NULL,
    conditionIdValue varchar(255) NOT NULL,
    conditionIdType varchar(255) NOT NULL,
    clinicalSignificance varchar(255) NOT NULL,
    clinicalSignificanceCitations varchar(255) NOT NULL,
    commentOnClinicalSignificance varchar(500) NOT NULL,
    dateLastEvaluated varchar(255) NOT NULL,
    url varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationESP (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInESP varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    minorAlleleFrequencyPercent varchar(255) NOT NULL,
    alleleFrequency varchar(255) NOT NULL,
    aaAlleleFrequency varchar(255) NOT NULL,
    eaAlleleFrequency varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationExAC (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInExAC varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    alleleFrequency varchar(255) NOT NULL,
    alleleFrequencyAFR varchar(255) NOT NULL,
    alleleFrequencyAMR varchar(255) NOT NULL,
    alleleFrequencyEAS varchar(255) NOT NULL,
    alleleFrequencyFIN varchar(255) NOT NULL,
    alleleFrequencyNFE varchar(255) NOT NULL,
    alleleFrequencyOTH varchar(255) NOT NULL,
    alleleFrequencySAS varchar(255) NOT NULL,
    alleleNumberAFR varchar(255) NOT NULL,
    alleleNumberAMR varchar(255) NOT NULL,
    alleleNumberEAS varchar(255) NOT NULL,
    alleleNumberFIN varchar(255) NOT NULL,
    alleleNumberNFE varchar(255) NOT NULL,
    alleleNumberOTH varchar(255) NOT NULL,
    alleleNumberSAS varchar(255) NOT NULL,
    homozygousCountAFR varchar(255) NOT NULL,
    homozygousCountAMR varchar(255) NOT NULL,
    homozygousCountEAS varchar(255) NOT NULL,
    homozygousCountFIN varchar(255) NOT NULL,
    homozygousCountNFE varchar(255) NOT NULL,
    homozygousCountOTH varchar(255) NOT NULL,
    homozygousCountSAS varchar(255) NOT NULL,
    alleleCountAFR varchar(255) NOT NULL,
    alleleCountAMR varchar(255) NOT NULL,
    alleleCountEAS varchar(255) NOT NULL,
    alleleCountFIN varchar(255) NOT NULL,
    alleleCountNFE varchar(255) NOT NULL,
    alleleCountOTH varchar(255) NOT NULL,
    alleleCountSAS varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationExLOVD (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInExLOVD varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    cooccurrenceLR varchar(255) NOT NULL,
    sumFamilyLR varchar(255) NOT NULL,
    segregationLR varchar(255) NOT NULL,
    posteriorProbability varchar(255) NOT NULL,
    missenseAnalysisPriorProbability varchar(255) NOT NULL,
    combinedPriorProbability varchar(255) NOT NULL,
    iarcClass varchar(255) NOT NULL,
    literatureSource varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);

CREATE TABLE brcaAnnotationLOVD (
    id int NOT NULL AUTO_INCREMENT,
    brcaId int NOT NULL,
    variantInLOVD varchar(255) NOT NULL,
    bxId varchar(255) NOT NULL,
    dbId varchar(255) NOT NULL,
    hgvsCDNA varchar(255) NOT NULL,
    hgvsProtein varchar(255) NOT NULL,
    rna varchar(255) NOT NULL,
    variantEffect varchar(255) NOT NULL,
    variantFrequency varchar(255) NOT NULL,
    variantHaplotype varchar(255) NOT NULL,
    geneticOrigin varchar(255) NOT NULL,
    functionalAnalysisTechnique varchar(255) NOT NULL,
    functionalAnalysisResult varchar(255) NOT NULL,
    submitters varchar(1000) NOT NULL,
    individuals varchar(255) NOT NULL,
    PRIMARY KEY (id),
    FOREIGN KEY (brcaId) REFERENCES brca(id)
);


