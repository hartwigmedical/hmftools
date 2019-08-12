package com.hartwig.hmftools.vicc;

import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.BRCA;
import com.hartwig.hmftools.vicc.datamodel.BRCApart1;
import com.hartwig.hmftools.vicc.datamodel.BRCApart2;
import com.hartwig.hmftools.vicc.datamodel.Civic;
import com.hartwig.hmftools.vicc.datamodel.CivicAvatars;
import com.hartwig.hmftools.vicc.datamodel.CivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.CivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.CivicDescription;
import com.hartwig.hmftools.vicc.datamodel.CivicDisease;
import com.hartwig.hmftools.vicc.datamodel.CivicDrugs;
import com.hartwig.hmftools.vicc.datamodel.CivicError;
import com.hartwig.hmftools.vicc.datamodel.CivicEvidenceItems;
import com.hartwig.hmftools.vicc.datamodel.CivicLastCommentedOn;
import com.hartwig.hmftools.vicc.datamodel.CivicLastModified;
import com.hartwig.hmftools.vicc.datamodel.CivicLastReviewed;
import com.hartwig.hmftools.vicc.datamodel.CivicLifecycleActions;
import com.hartwig.hmftools.vicc.datamodel.CivicOrganization;
import com.hartwig.hmftools.vicc.datamodel.CivicProfileImage;
import com.hartwig.hmftools.vicc.datamodel.CivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.CivicSource;
import com.hartwig.hmftools.vicc.datamodel.CivicUser;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.CivicVariantTypes;
import com.hartwig.hmftools.vicc.datamodel.CivicVariants;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivic;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicAvatars;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicClinicalTrial;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicCoordinates;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicDescription;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicDisease;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicDrugs;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicError;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicEvidenceItems;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicLastCommentedOn;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicLastModified;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicLastReviewed;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicLifecycleActions;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicOrganization;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicProfileImage;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicPublicationDate;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicSource;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicUser;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicVariantGroup;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicVariantTypes;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCivicVariants;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrials;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAreg;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchExonInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchExonsBoundries;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchFusions;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchParents;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchPositions;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTrialsTags;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchWGSaMap;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchWGSadataLocation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncokbGene;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncokbVariant;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePmkbGene;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePmkbVariant;
import com.hartwig.hmftools.vicc.datamodel.JaxTrials;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.JaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAreg;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstRightLeft;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstRightRight;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchBreg;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchExonInfo;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchExonsBoundries;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchFusionData;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchFusions;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchParents;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchPositions;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrials;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsContact;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsGeo;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsIntervation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsLocation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsOverallContact;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTrialsTags;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchVariantInfo;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchWGSaMap;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchWGSadataLocation;
import com.hartwig.hmftools.vicc.datamodel.OncoKbBiological;
import com.hartwig.hmftools.vicc.datamodel.Cgi;
import com.hartwig.hmftools.vicc.datamodel.OncoKbClinical;
import com.hartwig.hmftools.vicc.datamodel.OncoKbConsequence;
import com.hartwig.hmftools.vicc.datamodel.OncoKbDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.OncokbGene;
import com.hartwig.hmftools.vicc.datamodel.PmkbGene;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCA;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCApart1;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCApart2;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCgi;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJax;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxIndications;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxReferences;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTherapy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableOncokb;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePmkb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableSage;
import com.hartwig.hmftools.vicc.datamodel.ImmutableSequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTaxonomy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Jax;
import com.hartwig.hmftools.vicc.datamodel.JaxIndications;
import com.hartwig.hmftools.vicc.datamodel.JaxMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.JaxReferences;
import com.hartwig.hmftools.vicc.datamodel.JaxTherapy;
import com.hartwig.hmftools.vicc.datamodel.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.PmkbTissue;
import com.hartwig.hmftools.vicc.datamodel.Sage;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.PmkbTumor;
import com.hartwig.hmftools.vicc.datamodel.OncokbVariant;
import com.hartwig.hmftools.vicc.datamodel.PmkbVariant;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccJsonReader {
    private static final Logger LOGGER = LogManager.getLogger(ViccJsonReader.class);

    // SAGE records hold 8 field (no "feature names") while all other knowledgebases hold 9 records.
    private static final List<Integer> EXPECTED_VICC_ENTRY_SIZES = Lists.newArrayList(8, 9);

    private static final List<Integer> EXPECTED_ASSOCIATION_ELEMENT_SIZES = Lists.newArrayList(4, 5, 6, 7, 8, 9, 10, 11);
    private static final List<Integer> EXPECTED_FEATURES_ELEMENT_SIZES = Lists.newArrayList(2, 3, 5, 6, 8, 9, 10, 11, 12, 13, 14, 15, 16);
    private static final List<Integer> EXPECTED_SEQUENCE_ONTOLOGY_ELEMENT_SIZES = Lists.newArrayList(4, 5);
    private static final List<Integer> EXPECTED_GENE_IDENTIFIERS_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_EVIDENCE_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_EVIDENCE_INFO_ELEMENT_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES = Lists.newArrayList(1, 2);
    private static final List<Integer> EXPECTED_PHENOTYPE_ELEMENT_SIZES = Lists.newArrayList(2, 3, 4);
    private static final List<Integer> EXPECTED_PHENOTYPE_TYPE_ELEMENT_SIZES = Lists.newArrayList(3);

    private static final List<Integer> EXPECTED_CGI_ELEMENT_SIZES = Lists.newArrayList(23);

    private static final List<Integer> EXPECTED_BRCA_ELEMENT_SIZES = Lists.newArrayList(137);

    private static final List<Integer> EXPECTED_SAGE_ELEMENT_SIZES = Lists.newArrayList(8);

    private static final List<Integer> EXPECTED_PMKB_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_PMKB_TUMOR_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_PMKB_TISSUE_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_PMKB_VARIANT_ELEMENT_SIZES = Lists.newArrayList(21);
    private static final List<Integer> EXPECTED_PMKB_GENE_ELEMENT_SIZES = Lists.newArrayList(7);

    private static final List<Integer> EXPECTED_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES = Lists.newArrayList(11);
    private static final List<Integer> EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_ONCOKB_GENE_ELEMENT_SIZES = Lists.newArrayList(8);

    private static final List<Integer> EXPECTED_JAX_ELEMENT_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_JAX_THERAPY_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_INDICATIONS_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_JAX_REFERENCES_ELEMENT_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES = Lists.newArrayList(2);

    private static final List<Integer> EXPECTED_JAX_TRIALS_ELEMENT_SIZES = Lists.newArrayList(11);
    private static final List<Integer> EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES = Lists.newArrayList(2);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_ELEMENT_SIZES = Lists.newArrayList(34, 35, 36, 37, 38, 39, 40, 41, 42);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_AST_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES = Lists.newArrayList(3, 29, 30, 31);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES = Lists.newArrayList(8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES = Lists.newArrayList(3, 11);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES = Lists.newArrayList(13, 14, 16, 17, 18, 19);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PREVELANCE_SIZES = Lists.newArrayList(4, 6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_SOURCE_SIZES = Lists.newArrayList(8, 9, 10, 11, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TAGS_SIZES = Lists.newArrayList(3, 8, 9, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES = Lists.newArrayList(9, 14, 15, 16);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES = Lists.newArrayList(6);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES = Lists.newArrayList(10);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_PARENTS_SIZES = Lists.newArrayList(3, 4);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONS_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSADATA_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES = Lists.newArrayList(7, 9);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES =
            Lists.newArrayList(20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 39, 40, 41, 45);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES = Lists.newArrayList(3, 7);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_EXONSBOUNDRIES_SIZES =
            Lists.newArrayList(1, 2, 3, 5, 6, 7, 8, 9, 11, 13, 16, 17, 20, 21, 22, 24, 26, 27, 28, 29, 38, 41);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_POSITIONS_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES = Lists.newArrayList(1, 2, 15, 17);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_BREG_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_AREG_SIZES = Lists.newArrayList(2);

    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES = Lists.newArrayList(13, 14);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES = Lists.newArrayList(1, 2, 3, 4, 5);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES =
            Lists.newArrayList(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES =
            Lists.newArrayList(0, 1, 2, 3, 4, 5, 7, 8, 9, 10);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES = Lists.newArrayList(7, 8, 9, 10, 11, 12, 13);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES = Lists.newArrayList(0, 1, 2, 3);

    private static final List<Integer> EXPECTED_CIVIC_ELEMENT_SIZES = Lists.newArrayList(21);
    private static final List<Integer> EXPECTED_CIVIC_AVATARS_SIZES = Lists.newArrayList(4);
    private static final List<Integer> EXPECTED_CIVIC_COORDINATES_SIZES = Lists.newArrayList(12);
    private static final List<Integer> EXPECTED_CIVIC_DISEASES_SIZES = Lists.newArrayList(5);
    private static final List<Integer> EXPECTED_CIVIC_DRUGS_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES = Lists.newArrayList(17);
    private static final List<Integer> EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LAST_MODIFIED_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LAST_REVIEWED_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES = Lists.newArrayList(0, 3);
    private static final List<Integer> EXPECTED_CIVIC_ORGANIZATION_SIZES = Lists.newArrayList(0, 5);
    private static final List<Integer> EXPECTED_CIVIC_PROFILE_IMAGE_SIZES = Lists.newArrayList(5);
    private static final List<Integer> EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES = Lists.newArrayList(0, 1, 2, 3);
    private static final List<Integer> EXPECTED_CIVIC_SOURCE_SIZES = Lists.newArrayList(13);
    private static final List<Integer> EXPECTED_CIVIC_USER_SIZES = Lists.newArrayList(21);
    private static final List<Integer> EXPECTED_CIVIC_VARIANTTYPES_SIZES = Lists.newArrayList(6);
    private static final List<Integer> EXPECTED_CIVIC_VARIANTGROUPS_SIZES = Lists.newArrayList(5);
    private static final List<Integer> EXPECTED_CIVIC_VARIANTS_SIZES = Lists.newArrayList(10);
    private static final List<Integer> EXPECTED_CIVIC_PROVISIONALVALUES_SIZES = Lists.newArrayList(0, 1);
    private static final List<Integer> EXPECTED_CIVIC_DESCRIPTION_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_CIVIC_ERROR_SIZES = Lists.newArrayList(0);
    private static final List<Integer> EXPECTED_CIVIC_CLINICALTRIAL_SIZES = Lists.newArrayList(4);

    private ViccJsonReader() {
    }

    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);
        int number = 1;
        List<ViccEntry> entries = Lists.newArrayList();
        LOGGER.info("Reading VICC knowledgebase from " + jsonPath);
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryObject = parser.parse(reader).getAsJsonObject();
            //  LOGGER.info(number);
            number++;
            if (!EXPECTED_VICC_ENTRY_SIZES.contains(viccEntryObject.size())) {
                LOGGER.warn("Found " + viccEntryObject.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_VICC_ENTRY_SIZES);
                LOGGER.warn(viccEntryObject);
            }

            ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
            viccEntryBuilder.source(viccEntryObject.getAsJsonPrimitive("source").getAsString());
            viccEntryBuilder.genes(jsonArrayToStringList(viccEntryObject.getAsJsonArray("genes")));

            viccEntryBuilder.geneIdentifiers(createGeneIdentifiers(viccEntryObject));

            if (viccEntryObject.has("feature_names")) {
                JsonElement featureNames = viccEntryObject.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(createFeatures(viccEntryObject));

            JsonObject elementAssociation = viccEntryObject.getAsJsonObject("association");
            Set<String> keysAssociation = elementAssociation.getAsJsonObject().keySet();

            if (!EXPECTED_ASSOCIATION_ELEMENT_SIZES.contains(keysAssociation.size())) {
                LOGGER.warn("Found " + keysAssociation.size() + " in association rather than the expected "
                        + EXPECTED_ASSOCIATION_ELEMENT_SIZES);
                LOGGER.warn(keysAssociation);
            }

            viccEntryBuilder.association(createAssociation(elementAssociation));

            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("dev_tags")));

            JsonObject objectCgi = viccEntryObject.getAsJsonObject("cgi");
            if (viccEntryObject.has("cgi")) {
                Set<String> keysCgi = objectCgi.keySet();

                if (!EXPECTED_CGI_ELEMENT_SIZES.contains(keysCgi.size())) {
                    LOGGER.warn("Found " + keysCgi.size() + " in cgi rather than the expected " + EXPECTED_CGI_ELEMENT_SIZES);
                    LOGGER.warn(keysCgi);
                }
            }

            JsonObject objectBRCA = viccEntryObject.getAsJsonObject("brca");
            if (viccEntryObject.has("brca")) {
                Set<String> keysBRCA = objectBRCA.keySet();
                if (!EXPECTED_BRCA_ELEMENT_SIZES.contains(keysBRCA.size())) {
                    LOGGER.warn("Found " + keysBRCA.size() + " in brca rather than the expected " + EXPECTED_BRCA_ELEMENT_SIZES);
                    LOGGER.warn(keysBRCA);
                }
            }

            JsonObject objectSage = viccEntryObject.getAsJsonObject("sage");
            if (viccEntryObject.has("sage")) {
                Set<String> keysSage = objectSage.keySet();
                if (!EXPECTED_SAGE_ELEMENT_SIZES.contains(keysSage.size())) {
                    LOGGER.warn("Found " + keysSage.size() + " in sage rather than the expected " + EXPECTED_SAGE_ELEMENT_SIZES);
                    LOGGER.warn(keysSage);
                }
            }

            JsonObject objectPmkb = viccEntryObject.getAsJsonObject("pmkb");
            if (viccEntryObject.has("pmkb")) {
                Set<String> keysPmkb = objectPmkb.keySet();
                if (!EXPECTED_PMKB_ELEMENT_SIZES.contains(keysPmkb.size())) {
                    LOGGER.warn("Found " + keysPmkb.size() + " in pmkb rather than the expected " + EXPECTED_PMKB_ELEMENT_SIZES);
                    LOGGER.warn(keysPmkb);
                }
            }

            JsonObject objectOncokb = viccEntryObject.getAsJsonObject("oncokb");
            if (viccEntryObject.has("oncokb")) {
                Set<String> keysOncokb = objectOncokb.keySet();
                if (!EXPECTED_ONCOKB_ELEMENT_SIZES.contains(keysOncokb.size())) {
                    LOGGER.warn("Found " + keysOncokb.size() + " in oncokb rather than the expected " + EXPECTED_ONCOKB_ELEMENT_SIZES);
                    LOGGER.warn(keysOncokb);
                }
            }

            JsonObject objectJax = viccEntryObject.getAsJsonObject("jax");
            if (viccEntryObject.has("jax")) {
                Set<String> keysJax = objectJax.keySet();
                if (!EXPECTED_JAX_ELEMENT_SIZES.contains(keysJax.size())) {
                    LOGGER.warn("Found " + keysJax.size() + " in jax rather than the expected " + EXPECTED_JAX_ELEMENT_SIZES);
                    LOGGER.warn(keysJax);
                }
            }

            JsonObject objectJaxTrials = viccEntryObject.getAsJsonObject("jax_trials");

            if (viccEntryObject.has("jax_trials")) {
                Set<String> keysJaxTrials = objectJaxTrials.keySet();
                if (!EXPECTED_JAX_TRIALS_ELEMENT_SIZES.contains(keysJaxTrials.size())) {
                    LOGGER.warn("Found " + keysJaxTrials.size() + " in jax trials rather than the expected "
                            + EXPECTED_JAX_TRIALS_ELEMENT_SIZES);
                    LOGGER.warn(keysJaxTrials);
                }
            }

            JsonObject objectMolecularMatch = viccEntryObject.getAsJsonObject("molecularmatch");
            if (viccEntryObject.has("molecularmatch")) {
                Set<String> keysMolecularMatch = objectMolecularMatch.keySet();
                if (!EXPECTED_MOLECULARMATCH_ELEMENT_SIZES.contains(keysMolecularMatch.size())) {
                    LOGGER.warn("Found " + keysMolecularMatch.size() + " in molecular match rather than the expected "
                            + EXPECTED_MOLECULARMATCH_ELEMENT_SIZES);
                    LOGGER.warn(keysMolecularMatch);
                }
            }

            JsonObject objectMolecularMatchTrials = viccEntryObject.getAsJsonObject("molecularmatch_trials");
            if (viccEntryObject.has("molecularmatch_trials")) {
                Set<String> keysMolecularMatchTrials = objectMolecularMatchTrials.keySet();
                if (!EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES.contains(keysMolecularMatchTrials.size())) {
                    LOGGER.warn("Found " + keysMolecularMatchTrials.size() + " in molecular match trials rather than the expected "
                            + EXPECTED_MOLECULARMATCH_TRAILS_ELEMENT_SIZES);
                    LOGGER.warn(keysMolecularMatchTrials);
                }
            }

            JsonObject objectCivic = viccEntryObject.getAsJsonObject("civic");
            if (viccEntryObject.has("civic")) {
                Set<String> keysCivic = objectCivic.keySet();
                if (!EXPECTED_CIVIC_ELEMENT_SIZES.contains(keysCivic.size())) {
                    LOGGER.warn("Found " + keysCivic.size() + " in civic rather than the expected " + EXPECTED_CIVIC_ELEMENT_SIZES);
                    LOGGER.warn(keysCivic);
                }
            }

            if (viccEntryObject.has("cgi")) {
                viccEntryBuilder.KbSpecificObject(createCgi(objectCgi));
            } else if (viccEntryObject.has("brca")) {
                viccEntryBuilder.KbSpecificObject(createBRCA(objectBRCA));
            } else if (viccEntryObject.has("sage")) {
                viccEntryBuilder.KbSpecificObject(createSage(objectSage));
            } else if (viccEntryObject.has("pmkb")) {
                viccEntryBuilder.KbSpecificObject(createPmkb(objectPmkb));
            } else if (viccEntryObject.has("oncokb")) {
                if (objectOncokb.has("biological")) {
                    viccEntryBuilder.KbSpecificObject(createOncoKbBiological(objectOncokb));
                } else if (objectOncokb.has("clinical")) {
                    viccEntryBuilder.KbSpecificObject(createOncoKbClinical(objectOncokb));
                }
            } else if (viccEntryObject.has("jax")) {
                viccEntryBuilder.KbSpecificObject(createJax(objectJax));
            } else if (viccEntryObject.has("jax_trials")) {
                viccEntryBuilder.KbSpecificObject(createJaxTrials(objectJaxTrials));
            } else if (viccEntryObject.has("molecularmatch")) {
                viccEntryBuilder.KbSpecificObject(createMolecularMatch(objectMolecularMatch));
            } else if (viccEntryObject.has("molecularmatch_trials")) {
                viccEntryBuilder.KbSpecificObject(createMolecularMatchTrials(objectMolecularMatchTrials));
            } else if (viccEntryObject.has("civic")) {
                viccEntryBuilder.KbSpecificObject(createCivic(objectCivic));
            }
            entries.add(viccEntryBuilder.build());

        }
        reader.close();

        return entries;
    }

    @NotNull
    private static Civic createCivic(@NotNull JsonObject objectCivic) {
        return ImmutableCivic.builder()
                .variantGroups(!objectCivic.has("variant_groups")
                        ? null
                        : createVariantsGroups(objectCivic.getAsJsonArray("variant_groups")))
                .entrezName(objectCivic.getAsJsonPrimitive("entrez_name").getAsString())
                .variantTypes(createVariantTypes(objectCivic.getAsJsonArray("variant_types")))
                .civicActionabilityScore(objectCivic.get("civic_actionability_score").isJsonNull()
                        ? null
                        : objectCivic.getAsJsonPrimitive("civic_actionability_score").getAsString())
                .clinvarEntries(jsonArrayToStringList(objectCivic.getAsJsonArray("clinvar_entries")))
                .lifecycleActions(createLifeCycleActions(objectCivic.getAsJsonObject("lifecycle_actions")))
                .variantAliases(jsonArrayToStringList(objectCivic.getAsJsonArray("variant_aliases")))
                .alleleRegistryId(objectCivic.get("allele_registry_id").isJsonNull()
                        ? null
                        : objectCivic.getAsJsonPrimitive("allele_registry_id").getAsString())
                .provisional_values(createCivicDesription(objectCivic.getAsJsonObject("provisional_values")))
                .geneId(objectCivic.getAsJsonPrimitive("gene_id").getAsString())
                .name(objectCivic.getAsJsonPrimitive("name").getAsString())
                .evidenceItem(createEvidenceitems(objectCivic.getAsJsonArray("evidence_items")))
                .sources(createCivicSource(objectCivic.getAsJsonArray("sources")))
                .entrezId(objectCivic.getAsJsonPrimitive("entrez_id").getAsString())
                .assertions(jsonArrayToStringList(objectCivic.getAsJsonArray("assertions")))
                .hgvs_expressions(jsonArrayToStringList(objectCivic.getAsJsonArray("hgvs_expressions")))
                .errors(createCivicError(objectCivic.getAsJsonObject("errors")))
                .coordinates(createCoordinates(objectCivic.getAsJsonObject("coordinates")))
                .type(objectCivic.getAsJsonPrimitive("type").getAsString())
                .id(objectCivic.getAsJsonPrimitive("id").getAsString())
                .description(objectCivic.getAsJsonPrimitive("description").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicSource> createCivicSource(@NotNull JsonArray arraySource) {
        List<CivicSource> sourceList = Lists.newArrayList();
        for (JsonElement source : arraySource) {
            Set<String> keysSource = source.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_SOURCE_SIZES.contains(keysSource.size())) {
                LOGGER.warn("Found " + keysSource.size() + " in civic source rather than the expected " + EXPECTED_CIVIC_SOURCE_SIZES);
                LOGGER.warn(keysSource);
            }
            sourceList.add(ImmutableCivicSource.builder()
                    .status(source.getAsJsonObject().getAsJsonPrimitive("status").getAsString())
                    .openAccess(source.getAsJsonObject().get("open_access").isJsonNull()
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("open_access").getAsString())
                    .name(source.getAsJsonObject().get("name").isJsonNull()
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .journal(source.getAsJsonObject().get("journal").isJsonNull()
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("journal").getAsString())
                    .citation(source.getAsJsonObject().getAsJsonPrimitive("citation").getAsString())
                    .pmc_Id(source.getAsJsonObject().get("pmc_id").isJsonNull()
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("pmc_id").getAsString())
                    .fullJournalTitle(source.getAsJsonObject().get("full_journal_title").isJsonNull()
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("full_journal_title").getAsString())
                    .sourceUrl(source.getAsJsonObject().getAsJsonPrimitive("source_url").getAsString())
                    .clinicalTrials(createCivicClinicalTrials(source.getAsJsonObject().getAsJsonArray("clinical_trials")))
                    .pubmedId(source.getAsJsonObject().getAsJsonPrimitive("pubmed_id").getAsString())
                    .isReview(source.getAsJsonObject().getAsJsonPrimitive("is_review").getAsString())
                    .publicationDate(createPublicationDate(source.getAsJsonObject().getAsJsonObject("publication_date")))
                    .id(source.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .build());
        }

        return sourceList;
    }

    @NotNull
    private static CivicError createCivicError(@NotNull JsonObject objectError) {
        Set<String> keysError = objectError.getAsJsonObject().keySet();
        if (!EXPECTED_CIVIC_ERROR_SIZES.contains(keysError.size())) {
            LOGGER.warn("Found " + keysError.size() + " in civic error rather than the expected " + EXPECTED_CIVIC_ERROR_SIZES);
            LOGGER.warn(keysError);
        }
        return ImmutableCivicError.builder().build();
    }

    @NotNull
    private static CivicDescription createCivicDesription(@NotNull JsonObject objectProvisionalValues) {
        Set<String> keysProvisionalValues = objectProvisionalValues.getAsJsonObject().keySet();
        if (!EXPECTED_CIVIC_PROVISIONALVALUES_SIZES.contains(keysProvisionalValues.size())) {
            LOGGER.warn("Found " + keysProvisionalValues.size() + " in civic provisional values rather than the expected "
                    + EXPECTED_CIVIC_PROVISIONALVALUES_SIZES);
            LOGGER.warn(keysProvisionalValues);
        }

        JsonObject description =
                !objectProvisionalValues.has("description") ? null : objectProvisionalValues.getAsJsonObject("description");
        if (description != null) {
            Set<String> keysDescription = description.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_DESCRIPTION_SIZES.contains(keysDescription.size())) {
                LOGGER.warn("Found " + keysDescription.size() + " in civic description rather than the expected "
                        + EXPECTED_CIVIC_DESCRIPTION_SIZES);
                LOGGER.warn(keysDescription);
            }
            return ImmutableCivicDescription.builder()
                    .revision_id(!description.has("revision_id") || description.get("revision_id").isJsonNull()
                            ? null
                            : description.getAsJsonPrimitive("revision_id").getAsString())
                    .value(!description.has("value") || description.get("value").isJsonNull()
                            ? null
                            : description.getAsJsonPrimitive("value").getAsString())
                    .build();
        } else {
            return ImmutableCivicDescription.builder().revision_id("").value("").build();
        }

    }

    @NotNull
    private static List<CivicVariantGroup> createVariantsGroups(@NotNull JsonArray arrayVariantsGroup) {
        List<CivicVariantGroup> civicVariantGroupList = Lists.newArrayList();
        for (JsonElement variantGroup : arrayVariantsGroup) {
            Set<String> keysVariantGroup = variantGroup.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_VARIANTGROUPS_SIZES.contains(keysVariantGroup.size())) {
                LOGGER.warn("Found " + keysVariantGroup.size() + " in civic variant group rather than the expected "
                        + EXPECTED_CIVIC_VARIANTGROUPS_SIZES);
                LOGGER.warn(keysVariantGroup);
            }
            civicVariantGroupList.add(ImmutableCivicVariantGroup.builder()
                    .id(variantGroup.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .variants(createVariants(variantGroup.getAsJsonObject().getAsJsonArray("variants")))
                    .type(variantGroup.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .description(variantGroup.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .name(variantGroup.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return civicVariantGroupList;
    }

    @NotNull
    private static List<CivicVariants> createVariants(@NotNull JsonArray arrayVariants) {
        List<CivicVariants> variantsList = Lists.newArrayList();
        for (JsonElement variants : arrayVariants) {
            Set<String> keysVariants = variants.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_VARIANTS_SIZES.contains(keysVariants.size())) {
                LOGGER.warn(
                        "Found " + keysVariants.size() + " in civic variants rather than the expected " + EXPECTED_CIVIC_VARIANTS_SIZES);
                LOGGER.warn(keysVariants);
            }

            variantsList.add(ImmutableCivicVariants.builder()
                    .entrez_name(variants.getAsJsonObject().getAsJsonPrimitive("entrez_name").getAsString())
                    .variant_types(createVariantTypes(variants.getAsJsonObject().getAsJsonArray("variant_types")))
                    .description(variants.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .civic_actionability_score(!variants.getAsJsonObject().has("civic_actionability_score") || variants.getAsJsonObject()
                            .get("civic_actionability_score")
                            .isJsonNull() ? null : variants.getAsJsonObject().getAsJsonPrimitive("civic_actionability_score").getAsString())
                    .gene_id(variants.getAsJsonObject().getAsJsonPrimitive("gene_id").getAsString())
                    .entrez_id(variants.getAsJsonObject().getAsJsonPrimitive("entrez_id").getAsString())
                    .coordinates(!variants.getAsJsonObject().has("createCoordinates")
                            ? null
                            : createCoordinates(variants.getAsJsonObject().getAsJsonObject("createCoordinates")))
                    .type(variants.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(variants.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(variants.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return variantsList;
    }

    @NotNull
    private static CivicCoordinates createCoordinates(@NotNull JsonObject objectCoordinates) {
        Set<String> keysCoordinates = objectCoordinates.keySet();
        if (!EXPECTED_CIVIC_COORDINATES_SIZES.contains(keysCoordinates.size())) {
            LOGGER.warn("Found " + keysCoordinates.size() + " in civic coordinates rather than the expectedd "
                    + EXPECTED_CIVIC_COORDINATES_SIZES);
            LOGGER.warn(keysCoordinates);
        }

        return ImmutableCivicCoordinates.builder()
                .chromosome2(objectCoordinates.get("chromosome2").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("chromosome2").getAsString())
                .referenceBases(objectCoordinates.get("reference_bases").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("reference_bases").getAsString())
                .start2(objectCoordinates.get("start2").isJsonNull() ? null : objectCoordinates.getAsJsonPrimitive("start2").getAsString())
                .variantBases(objectCoordinates.get("variant_bases").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("variant_bases").getAsString())
                .stop(objectCoordinates.get("stop").isJsonNull() ? null : objectCoordinates.getAsJsonPrimitive("stop").getAsString())
                .stop2(objectCoordinates.get("stop2").isJsonNull() ? null : objectCoordinates.getAsJsonPrimitive("stop2").getAsString())
                .representativeTranscript2(objectCoordinates.get("representative_transcript2").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("representative_transcript2").getAsString())
                .start(objectCoordinates.get("start").isJsonNull() ? null : objectCoordinates.getAsJsonPrimitive("start").getAsString())
                .representativeTranscript(objectCoordinates.get("representative_transcript").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("representative_transcript").getAsString())
                .ensemblVersion(objectCoordinates.get("ensembl_version").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("ensembl_version").getAsString())
                .chromosome(objectCoordinates.get("chromosome").isJsonNull()
                        ? null
                        : objectCoordinates.getAsJsonPrimitive("chromosome").getAsString())
                .referenceBuild("")
                .build();
    }

    @NotNull
    private static List<CivicEvidenceItems> createEvidenceitems(@NotNull JsonArray evidenceItemsArray) {
        List<CivicEvidenceItems> evidenceItemsList = Lists.newArrayList();
        for (JsonElement evideneItem : evidenceItemsArray) {
            Set<String> keysEvidenceItems = evideneItem.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES.contains(keysEvidenceItems.size())) {
                LOGGER.warn("Found " + keysEvidenceItems.size() + " in civic evidence items rather than the expected "
                        + EXPECTED_CIVIC_EVIDENCE_ITEMS_SIZES);
                LOGGER.warn(keysEvidenceItems);
            }

            evidenceItemsList.add(ImmutableCivicEvidenceItems.builder()
                    .status(evideneItem.getAsJsonObject().getAsJsonPrimitive("status").getAsString())
                    .rating(evideneItem.getAsJsonObject().get("rating").isJsonNull()
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("rating").getAsString())
                    .drugInteractionType(evideneItem.getAsJsonObject().get("drug_interaction_type").isJsonNull()
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("drug_interaction_type").getAsString())
                    .description(evideneItem.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .openChangeCount(evideneItem.getAsJsonObject().getAsJsonPrimitive("open_change_count").getAsString())
                    .evidenceType(evideneItem.getAsJsonObject().getAsJsonPrimitive("evidence_type").getAsString())
                    .drugs(createDrugs(evideneItem.getAsJsonObject().getAsJsonArray("drugs")))
                    .variantOrigin(evideneItem.getAsJsonObject().get("variant_origin").isJsonNull()
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("variant_origin").getAsString())
                    .disease(createDiseases(evideneItem.getAsJsonObject().getAsJsonObject("disease")))
                    .source(createSource(evideneItem.getAsJsonObject().getAsJsonObject("source")))
                    .evidenceDirection(evideneItem.getAsJsonObject().get("evidence_direction").isJsonNull()
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("evidence_direction").toString())
                    .variantId(!evideneItem.getAsJsonObject().has("variant_id")
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("variant_id").getAsString())
                    .clinicalSignificance(evideneItem.getAsJsonObject().get("clinical_significance").isJsonNull()
                            ? null
                            : evideneItem.getAsJsonObject().getAsJsonPrimitive("clinical_significance").getAsString())
                    .evidenceLevel(evideneItem.getAsJsonObject().getAsJsonPrimitive("evidence_level").getAsString())
                    .type(evideneItem.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(evideneItem.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(evideneItem.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return evidenceItemsList;
    }

    @NotNull
    private static CivicSource createSource(@NotNull JsonObject objectSource) {
        Set<String> keysSource = objectSource.keySet();
        if (!EXPECTED_CIVIC_SOURCE_SIZES.contains(keysSource.size())) {
            LOGGER.warn("Found " + keysSource.size() + " in civic source rather than the expected " + EXPECTED_CIVIC_SOURCE_SIZES);
            LOGGER.warn(keysSource);
        }
        return ImmutableCivicSource.builder()
                .status(objectSource.getAsJsonPrimitive("status").getAsString())
                .openAccess(objectSource.get("open_access").isJsonNull()
                        ? null
                        : objectSource.getAsJsonPrimitive("open_access").getAsString())
                .name(objectSource.get("name").isJsonNull() ? null : objectSource.getAsJsonPrimitive("name").getAsString())
                .journal(objectSource.get("journal").isJsonNull() ? null : objectSource.getAsJsonPrimitive("journal").getAsString())
                .citation(objectSource.getAsJsonPrimitive("citation").getAsString())
                .pmc_Id(objectSource.get("pmc_id").isJsonNull() ? null : objectSource.getAsJsonPrimitive("pmc_id").getAsString())
                .fullJournalTitle(objectSource.get("full_journal_title").isJsonNull()
                        ? null
                        : objectSource.getAsJsonPrimitive("full_journal_title").getAsString())
                .sourceUrl(objectSource.getAsJsonPrimitive("source_url").getAsString())
                .clinicalTrials(createCivicClinicalTrials(objectSource.getAsJsonArray("clinical_trials")))
                .pubmedId(objectSource.getAsJsonPrimitive("pubmed_id").getAsString())
                .isReview(objectSource.getAsJsonPrimitive("is_review").getAsString())
                .publicationDate(createPublicationDate(objectSource.getAsJsonObject("publication_date")))
                .id(objectSource.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicClinicalTrial> createCivicClinicalTrials(@NotNull JsonArray arrayClinicalTrial) {
        List<CivicClinicalTrial> clinicalTrialList = Lists.newArrayList();
        for (JsonElement clinicalTrial : arrayClinicalTrial) {
            Set<String> keysClinicalTrials = clinicalTrial.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_CLINICALTRIAL_SIZES.contains(keysClinicalTrials.size())) {
                LOGGER.warn("Found " + keysClinicalTrials.size() + " in civic clinial trial rather than the expected "
                        + EXPECTED_CIVIC_CLINICALTRIAL_SIZES);
                LOGGER.warn(keysClinicalTrials);
            }

            clinicalTrialList.add(ImmutableCivicClinicalTrial.builder()
                    .nct_id(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("nct_id").getAsString())
                    .description(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .clinical_trial_url(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("clinical_trial_url").getAsString())
                    .name(clinicalTrial.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return clinicalTrialList;
    }

    @NotNull
    private static CivicPublicationDate createPublicationDate(@NotNull JsonObject objectPublicationDate) {
        Set<String> keysPublicationDate = objectPublicationDate.keySet();
        if (!EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES.contains(keysPublicationDate.size())) {
            LOGGER.warn("Found " + keysPublicationDate.size() + " in civic publication date rather than the expected "
                    + EXPECTED_CIVIC_PUBLICATIONS_DATE_SIZES);
            LOGGER.warn(keysPublicationDate);
        }

        return ImmutableCivicPublicationDate.builder()
                .year(!objectPublicationDate.has("year") ? null : objectPublicationDate.getAsJsonPrimitive("year").getAsString())
                .day(!objectPublicationDate.has("day") ? null : objectPublicationDate.getAsJsonPrimitive("day").getAsString())
                .month(!objectPublicationDate.has("month") ? null : objectPublicationDate.getAsJsonPrimitive("month").getAsString())
                .build();
    }

    @NotNull
    private static CivicDisease createDiseases(@NotNull JsonObject objectDisease) {
        Set<String> keysDisease = objectDisease.keySet();
        if (!EXPECTED_CIVIC_DISEASES_SIZES.contains(keysDisease.size())) {
            LOGGER.warn("Found " + keysDisease.size() + " in civic disease rather than the expected " + EXPECTED_CIVIC_DISEASES_SIZES);
            LOGGER.warn(keysDisease);
        }

        return ImmutableCivicDisease.builder()
                .doid(objectDisease.get("doid").isJsonNull() ? null : objectDisease.getAsJsonPrimitive("doid").getAsString())
                .url(objectDisease.getAsJsonPrimitive("url").getAsString())
                .displayName(objectDisease.getAsJsonPrimitive("display_name").getAsString())
                .id(objectDisease.getAsJsonPrimitive("id").getAsString())
                .name(objectDisease.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicDrugs> createDrugs(@NotNull JsonArray arrayDrugs) {
        List<CivicDrugs> drugsList = Lists.newArrayList();
        for (JsonElement drug : arrayDrugs) {
            Set<String> keysDrugs = drug.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_DRUGS_SIZES.contains(keysDrugs.size())) {
                LOGGER.warn("Found " + keysDrugs.size() + " in civic drugs rather than the expected " + EXPECTED_CIVIC_DRUGS_SIZES);
                LOGGER.warn(keysDrugs);
            }

            drugsList.add(ImmutableCivicDrugs.builder()
                    .pubchemId(drug.getAsJsonObject().get("pubchem_id").isJsonNull()
                            ? null
                            : drug.getAsJsonObject().getAsJsonPrimitive("pubchem_id").getAsString())
                    .id(drug.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(drug.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return drugsList;
    }

    @NotNull
    private static CivicLifecycleActions createLifeCycleActions(@NotNull JsonObject objectLifeCycleActions) {
        Set<String> keysLifecycleActions = objectLifeCycleActions.keySet();
        if (!EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES.contains(keysLifecycleActions.size())) {
            LOGGER.warn("Found " + keysLifecycleActions.size() + " in civic lifecycle actions rather than the expected "
                    + EXPECTED_CIVIC_LIFECYCLE_ACTIONS_SIZES);
            LOGGER.warn(keysLifecycleActions);
        }

        return ImmutableCivicLifecycleActions.builder()
                .lastCommentedOn(objectLifeCycleActions.getAsJsonObject("last_commented_on") == null
                        ? null
                        : createLastCommentOn(objectLifeCycleActions.getAsJsonObject("last_commented_on")))
                .lastModified(objectLifeCycleActions.getAsJsonObject("last_modified") == null
                        ? null
                        : createLastModified(objectLifeCycleActions.getAsJsonObject("last_modified")))
                .lastReviewed(objectLifeCycleActions.getAsJsonObject("last_reviewed") == null
                        ? null
                        : createLastReviewed(objectLifeCycleActions.getAsJsonObject("last_reviewed")))
                .build();
    }

    @NotNull
    private static CivicLastCommentedOn createLastCommentOn(@NotNull JsonObject objectLastCommned) {
        Set<String> keysLastCommentedOn = objectLastCommned.keySet();
        if (!EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES.contains(keysLastCommentedOn.size())) {
            LOGGER.warn("Found " + keysLastCommentedOn.size() + " in civic last commented on rather than the expected "
                    + EXPECTED_CIVIC_LAST_COMMENTED_ON_SIZES);
            LOGGER.warn(keysLastCommentedOn);
        }

        return ImmutableCivicLastCommentedOn.builder()
                .timestamp(objectLastCommned.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastCommned.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicLastModified createLastModified(@NotNull JsonObject objectLastModified) {
        Set<String> keysLastModified = objectLastModified.keySet();
        if (!EXPECTED_CIVIC_LAST_MODIFIED_SIZES.contains(keysLastModified.size())) {
            LOGGER.warn("Found " + keysLastModified.size() + " in civic last modified rather than the expected"
                    + EXPECTED_CIVIC_LAST_MODIFIED_SIZES);
            LOGGER.warn(keysLastModified);
        }

        return ImmutableCivicLastModified.builder()
                .timestamp(objectLastModified.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastModified.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicLastReviewed createLastReviewed(@NotNull JsonObject objectLastReviewed) {
        Set<String> keysLastReviewed = objectLastReviewed.keySet();
        if (!EXPECTED_CIVIC_LAST_REVIEWED_SIZES.contains(keysLastReviewed.size())) {
            LOGGER.warn("Found " + keysLastReviewed.size() + " in civic last reviewed rather than the expected "
                    + EXPECTED_CIVIC_LAST_REVIEWED_SIZES);
            LOGGER.warn(keysLastReviewed);
        }

        return ImmutableCivicLastReviewed.builder()
                .timestamp(objectLastReviewed.getAsJsonPrimitive("timestamp").getAsString())
                .user(createCivicUser(objectLastReviewed.getAsJsonObject("user")))
                .build();
    }

    @NotNull
    private static CivicUser createCivicUser(@NotNull JsonObject objectUser) {
        Set<String> keysUser = objectUser.keySet();
        if (!EXPECTED_CIVIC_USER_SIZES.contains(keysUser.size())) {
            LOGGER.warn("Found " + keysUser.size() + " in civic user rather than the expected " + EXPECTED_CIVIC_USER_SIZES);
            LOGGER.warn(keysUser);
        }

        return ImmutableCivicUser.builder()
                .username(objectUser.getAsJsonPrimitive("username").getAsString())
                .areaOfExpertise(objectUser.get("area_of_expertise").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("area_of_expertise").getAsString())
                .organization(createOrganization(objectUser.getAsJsonObject("organization")))
                .twitterHandle(objectUser.get("twitter_handle").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("twitter_handle").getAsString())
                .name(objectUser.getAsJsonPrimitive("name").getAsString())
                .bio(objectUser.get("bio").isJsonNull() ? null : objectUser.getAsJsonPrimitive("bio").getAsString())
                .url(objectUser.get("url").isJsonNull() ? null : objectUser.getAsJsonPrimitive("url").getAsString())
                .createdAt(objectUser.getAsJsonPrimitive("created_at").getAsString())
                .avatars(createAvatars(objectUser.getAsJsonObject("avatars")))
                .acceptedLicense(objectUser.get("accepted_license").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("accepted_license").getAsString())
                .affiliation(objectUser.get("affiliation").isJsonNull() ? null : objectUser.getAsJsonPrimitive("affiliation").getAsString())
                .avatarUrl(objectUser.getAsJsonPrimitive("avatar_url").getAsString())
                .role(objectUser.getAsJsonPrimitive("role").getAsString())
                .facebookProfile(objectUser.get("facebook_profile").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("facebook_profile").getAsString())
                .linkedinProfile(objectUser.get("linkedin_profile").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("linkedin_profile").getAsString())
                .orcid(objectUser.get("orcid").isJsonNull() ? null : objectUser.getAsJsonPrimitive("orcid").getAsString())
                .displayName(objectUser.getAsJsonPrimitive("display_name").getAsString())
                .lastSeenAt(objectUser.get("last_seen_at").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("last_seen_at").getAsString())
                .featuredExpert(objectUser.getAsJsonPrimitive("featured_expert").getAsString())
                .id(objectUser.getAsJsonPrimitive("id").getAsString())
                .signupComplete(objectUser.get("signup_complete").isJsonNull()
                        ? null
                        : objectUser.getAsJsonPrimitive("signup_complete").getAsString())
                .build();
    }

    @NotNull
    private static CivicOrganization createOrganization(@NotNull JsonObject objectOrganization) {
        Set<String> keysOrganization = objectOrganization.keySet();
        if (!EXPECTED_CIVIC_ORGANIZATION_SIZES.contains(keysOrganization.size())) {
            LOGGER.warn("Found " + keysOrganization.size() + " in civic organization rather than the expected "
                    + EXPECTED_CIVIC_ORGANIZATION_SIZES);
            LOGGER.warn(keysOrganization);
        }

        return ImmutableCivicOrganization.builder()
                .url(!objectOrganization.has("url") ? null : objectOrganization.getAsJsonPrimitive("url").getAsString())
                .id(!objectOrganization.has("id") ? null : objectOrganization.getAsJsonPrimitive("id").getAsString())
                .profileImage(!objectOrganization.has("profile_image")
                        ? null
                        : createProfileImage(objectOrganization.getAsJsonObject("profile_image")))
                .description(!objectOrganization.has("description")
                        ? null
                        : objectOrganization.getAsJsonPrimitive("description").getAsString())
                .name(!objectOrganization.has("name") ? null : objectOrganization.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static CivicProfileImage createProfileImage(@NotNull JsonObject objectProfileImage) {
        Set<String> keysProfileImage = objectProfileImage.keySet();
        if (!EXPECTED_CIVIC_PROFILE_IMAGE_SIZES.contains(keysProfileImage.size())) {
            LOGGER.warn("Found " + keysProfileImage.size() + " in civic profile image rather than the expected "
                    + EXPECTED_CIVIC_PROFILE_IMAGE_SIZES);
            LOGGER.warn(keysProfileImage);
        }

        return ImmutableCivicProfileImage.builder()
                .x32(objectProfileImage.getAsJsonPrimitive("x32").getAsString())
                .x256(objectProfileImage.getAsJsonPrimitive("x256").getAsString())
                .x14(objectProfileImage.getAsJsonPrimitive("x14").getAsString())
                .x64(objectProfileImage.getAsJsonPrimitive("x64").getAsString())
                .x128(objectProfileImage.getAsJsonPrimitive("x128").getAsString())
                .build();
    }

    @NotNull
    private static CivicAvatars createAvatars(@NotNull JsonObject objectAvatars) {
        Set<String> keysAvatars = objectAvatars.keySet();
        if (!EXPECTED_CIVIC_AVATARS_SIZES.contains(keysAvatars.size())) {
            LOGGER.warn("Found " + keysAvatars.size() + " in civic avatars rather than the expected " + EXPECTED_CIVIC_AVATARS_SIZES);
            LOGGER.warn(keysAvatars);
        }

        return ImmutableCivicAvatars.builder()
                .x32(objectAvatars.getAsJsonPrimitive("x32").getAsString())
                .x14(objectAvatars.getAsJsonPrimitive("x14").getAsString())
                .x64(objectAvatars.getAsJsonPrimitive("x64").getAsString())
                .x128(objectAvatars.getAsJsonPrimitive("x128").getAsString())
                .build();
    }

    @NotNull
    private static List<CivicVariantTypes> createVariantTypes(@NotNull JsonArray arrayvariantTypes) {
        List<CivicVariantTypes> civicVariantTypesList = Lists.newArrayList();
        for (JsonElement variantTypes : arrayvariantTypes) {
            Set<String> keysVariantTypes = variantTypes.getAsJsonObject().keySet();
            if (!EXPECTED_CIVIC_VARIANTTYPES_SIZES.contains(keysVariantTypes.size())) {
                LOGGER.warn("Found " + keysVariantTypes.size() + " in civic variant types rather than the expected "
                        + EXPECTED_CIVIC_VARIANTTYPES_SIZES);
                LOGGER.warn(keysVariantTypes);
            }

            civicVariantTypesList.add(ImmutableCivicVariantTypes.builder()
                    .displayName(variantTypes.getAsJsonObject().getAsJsonPrimitive("display_name").getAsString())
                    .description(variantTypes.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .url(variantTypes.getAsJsonObject().getAsJsonPrimitive("url").getAsString())
                    .soId(variantTypes.getAsJsonObject().getAsJsonPrimitive("so_id").getAsString())
                    .id(variantTypes.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(variantTypes.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return civicVariantTypesList;
    }

    @NotNull
    private static MolecularMatchTrials createMolecularMatchTrials(@NotNull JsonObject objectMolecularMatchTrials) {
        return ImmutableMolecularMatchTrials.builder()
                .status(objectMolecularMatchTrials.getAsJsonPrimitive("status").getAsString())
                .startDate(!objectMolecularMatchTrials.has("startDate") || objectMolecularMatchTrials.get("startDate").isJsonNull()
                        ? null
                        : objectMolecularMatchTrials.getAsJsonPrimitive("startDate").getAsString())
                .title(objectMolecularMatchTrials.getAsJsonPrimitive("title").getAsString())
                .molecularAlterations(jsonArrayToStringList(objectMolecularMatchTrials.getAsJsonArray("molecularAlterations")))
                .score(objectMolecularMatchTrials.getAsJsonPrimitive("_score").getAsString())
                .intervation(createMolecularMatchTrialsIntervations(objectMolecularMatchTrials.getAsJsonArray("interventions")))
                .locations(createMolecularMatchTrialsLocations(objectMolecularMatchTrials.getAsJsonArray("locations")))
                .briefTitle(objectMolecularMatchTrials.get("briefTitle").isJsonNull()
                        ? null
                        : objectMolecularMatchTrials.getAsJsonPrimitive("briefTitle").getAsString())
                .overallContact(objectMolecularMatchTrials.get("overallContact").isJsonNull()
                        ? null
                        : createMolecularMatchTrialsOverallContact(objectMolecularMatchTrials.getAsJsonObject("overallContact")))
                .link(objectMolecularMatchTrials.getAsJsonPrimitive("link").getAsString())
                .phase(objectMolecularMatchTrials.getAsJsonPrimitive("phase").getAsString())
                .tags(createMolecularMatchTrialsTags(objectMolecularMatchTrials.getAsJsonArray("tags")))
                .id(objectMolecularMatchTrials.getAsJsonPrimitive("id").getAsString())
                .studyType(objectMolecularMatchTrials.getAsJsonPrimitive("studyType").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsLocations> createMolecularMatchTrialsLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchTrialsLocations> locationsList = Lists.newArrayList();
        for (JsonElement location : arrayLocations) {
            Set<String> keysLocation = location.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES.contains(keysLocation.size())) {
                LOGGER.warn("Found " + keysLocation.size() + " in molecular match trials locations rather than the expected"
                        + EXPECTED_MOLECULARMATCH_TRAILS_LOCATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysLocation);
            }

            locationsList.add(ImmutableMolecularMatchTrialsLocations.builder()
                    .status(location.getAsJsonObject().getAsJsonPrimitive("status").getAsString())
                    .last_name(!location.getAsJsonObject().has("last_name")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("last_name").getAsString())
                    .email(!location.getAsJsonObject().has("email")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("email").getAsString())
                    .phone(!location.getAsJsonObject().has("phone")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone").getAsString())
                    .phone_backup(!location.getAsJsonObject().has("phone_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_backup").getAsString())
                    .email_backup(!location.getAsJsonObject().has("email_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("email_backup").getAsString())
                    .last_name_backup(!location.getAsJsonObject().has("last_name_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("last_name_backup").getAsString())
                    .phone_ext_backup(!location.getAsJsonObject().has("phone_ext_backup")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_ext_backup").getAsString())
                    .phone_ext(!location.getAsJsonObject().has("phone_ext")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("phone_ext").getAsString())
                    .city(!location.getAsJsonObject().has("city")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("city").getAsString())
                    .valid(!location.getAsJsonObject().has("_valid")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("_valid").getAsString())
                    .zip(!location.getAsJsonObject().has("zip") ? null : location.getAsJsonObject().getAsJsonPrimitive("zip").getAsString())
                    .created(!location.getAsJsonObject().has("created")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("created").getAsString())
                    .country(!location.getAsJsonObject().has("country")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("country").getAsString())
                    .number(!location.getAsJsonObject().has("number")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("number").getAsString())
                    .id(!location.getAsJsonObject().has("id") ? null : location.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .lastUpdated(!location.getAsJsonObject().has("lastUpdated")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("lastUpdated").getAsString())
                    .contact(!location.getAsJsonObject().has("contact")
                            ? null
                            : createMolecularMatchTrialsContact(location.getAsJsonObject().getAsJsonObject("contact")))
                    .state(!location.getAsJsonObject().has("state")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("state").getAsString())
                    .street(!location.getAsJsonObject().has("street")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("street").getAsString())
                    .location(!location.getAsJsonObject().has("location") || location.getAsJsonObject().get("location").isJsonNull()
                            ? null
                            : createMolecularMatchTrialsLocation(location.getAsJsonObject().getAsJsonObject("location")))
                    .po_box(!location.getAsJsonObject().has("po_box")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("po_box").getAsString())
                    .failedGeocode(!location.getAsJsonObject().has("failedGeocode")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("failedGeocode").getAsString())
                    .geo(!location.getAsJsonObject().has("geo")
                            ? null
                            : createMolecularMatchTrialsGeo(location.getAsJsonObject().getAsJsonObject("geo")))
                    .validMessage(!location.getAsJsonObject().has("_validMessage")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("_validMessage").getAsString())
                    .name(!location.getAsJsonObject().has("name")
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<MolecularMatchTrialsIntervation> createMolecularMatchTrialsIntervations(@NotNull JsonArray intervationsArray) {
        List<MolecularMatchTrialsIntervation> molecularMatchTrialsIntervationList = Lists.newArrayList();
        for (JsonElement intervation : intervationsArray) {
            Set<String> keysClinical = intervation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES.contains(keysClinical.size())) {
                LOGGER.warn("Found " + keysClinical.size() + " in molecular match trials intervation rather than the expected "
                        + EXPECTED_MOLECULARMATCH_TRAILS_INTERVATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysClinical);
            }

            molecularMatchTrialsIntervationList.add(ImmutableMolecularMatchTrialsIntervation.builder()
                    .intervention_name(!intervation.getAsJsonObject().has("intervention_name")
                            ? null
                            : intervation.getAsJsonObject().getAsJsonPrimitive("intervention_name").getAsString())
                    .other_name(!intervation.getAsJsonObject().has("other_name")
                            ? null
                            : otherNameMolecularMatchTrials(intervation.getAsJsonObject()))
                    .description(!intervation.getAsJsonObject().has("description")
                            ? null
                            : intervation.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .arm_group_label(!intervation.getAsJsonObject().has("arm_group_label")
                            ? null
                            : armGroupLabelMolecularMatchTrials(intervation.getAsJsonObject()))
                    .intervention_type(!intervation.getAsJsonObject().has("intervention_type")
                            ? null
                            : intervation.getAsJsonObject().getAsJsonPrimitive("intervention_type").getAsString())
                    .build());
        }
        return molecularMatchTrialsIntervationList;
    }

    @NotNull
    private static Iterable<String> otherNameMolecularMatchTrials(@NotNull JsonObject otherNameObject) {
        if (otherNameObject.get("other_name").isJsonPrimitive()) {
            return Arrays.asList(otherNameObject.getAsJsonPrimitive("other_name").getAsString());
        } else {
            return jsonArrayToStringList(otherNameObject.getAsJsonArray("other_name"));
        }

    }

    @NotNull
    private static Iterable<String> armGroupLabelMolecularMatchTrials(@NotNull JsonObject armGroupLabel) {
        if (armGroupLabel.get("arm_group_label").isJsonArray()) {
            return jsonArrayToStringList(armGroupLabel.getAsJsonArray("arm_group_label"));
        } else if (armGroupLabel.get("arm_group_label").isJsonPrimitive()) {
            return Arrays.asList(armGroupLabel.getAsJsonPrimitive("arm_group_label").getAsString());
        } else {
            return Lists.newArrayList();
        }
    }

    @NotNull
    private static MolecularMatchTrialsOverallContact createMolecularMatchTrialsOverallContact(@NotNull JsonObject overallContactObject) {
        Set<String> keysOverallContact = overallContactObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES.contains(keysOverallContact.size())) {
            LOGGER.warn("Found " + keysOverallContact.size() + " in molecular match trials overall contact rather than the expected "
                    + EXPECTED_MOLECULARMATCH_TRAILS_OVERALL_CONTACT_ELEMENT_SIZES);
            LOGGER.warn(keysOverallContact);
        }
        return ImmutableMolecularMatchTrialsOverallContact.builder()
                .phone(!overallContactObject.has("phone") || overallContactObject.get("phone").isJsonNull()
                        ? null
                        : overallContactObject.getAsJsonPrimitive("phone").getAsString())
                .last_name(!overallContactObject.has("last_name")
                        ? null
                        : overallContactObject.getAsJsonPrimitive("last_name").getAsString())
                .email(!overallContactObject.has("email") ? null : overallContactObject.getAsJsonPrimitive("email").getAsString())
                .affiliation(!overallContactObject.has("affiliation") || overallContactObject.get("affiliation").isJsonNull()
                        ? null
                        : overallContactObject.getAsJsonPrimitive("affiliation").getAsString())
                .phone_ext(!overallContactObject.has("phone_ext")
                        ? null
                        : overallContactObject.getAsJsonPrimitive("phone_ext").getAsString())
                .country(!overallContactObject.has("country") ? null : overallContactObject.getAsJsonPrimitive("country").getAsString())
                .city(!overallContactObject.has("city") ? null : overallContactObject.getAsJsonPrimitive("city").getAsString())
                .name(!overallContactObject.has("name") ? null : overallContactObject.getAsJsonPrimitive("name").getAsString())
                .zip(!overallContactObject.has("zip") ? null : overallContactObject.getAsJsonPrimitive("zip").getAsString())
                .url(!overallContactObject.has("url") ? null : overallContactObject.getAsJsonPrimitive("url").getAsString())
                .street(!overallContactObject.has("street") ? null : overallContactObject.getAsJsonPrimitive("street").getAsString())
                .type(!overallContactObject.has("type") ? null : overallContactObject.getAsJsonPrimitive("type").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsGeo createMolecularMatchTrialsGeo(@NotNull JsonObject geoObject) {
        Set<String> keySGeo = geoObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES.contains(keySGeo.size())) {
            LOGGER.warn("Found " + keySGeo.size() + " in molecular match trials geo rather than the expected "
                    + EXPECTED_MOLECULARMATCH_TRAILS_GEO_ELEMENT_SIZES);
            LOGGER.warn(keySGeo);
        }
        return ImmutableMolecularMatchTrialsGeo.builder()
                .lat(geoObject.getAsJsonPrimitive("lat").getAsString())
                .lon(geoObject.getAsJsonPrimitive("lon").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsLocation createMolecularMatchTrialsLocation(@NotNull JsonObject locationObject) {
        Set<String> keysLocation = locationObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES.contains(keysLocation.size())) {
            LOGGER.warn("Found " + keysLocation.size() + " in molecular match trials location rather than the expected "
                    + EXPECTED_MOLECULARMATCH_TRAILS_LOCATION_ELEMENT_SIZES);
            LOGGER.warn(keysLocation);
        }
        return ImmutableMolecularMatchTrialsLocation.builder()
                .type(locationObject.getAsJsonPrimitive("type").getAsString())
                .coordinates(jsonArrayToStringList(locationObject.getAsJsonArray("coordinates")))
                .build();
    }

    @NotNull
    private static MolecularMatchTrialsContact createMolecularMatchTrialsContact(@NotNull JsonObject contactObject) {
        Set<String> keysContact = contactObject.getAsJsonObject().keySet();
        if (!EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES.contains(keysContact.size())) {
            LOGGER.warn("Found " + keysContact.size() + " in molecular match trials contact rather than the expected "
                    + EXPECTED_MOLECULARMATCH_TRAILS_CONTACT_ELEMENT_SIZES);
            LOGGER.warn(keysContact);
        }
        return ImmutableMolecularMatchTrialsContact.builder()
                .phone(!contactObject.has("phone") ? null : contactObject.getAsJsonPrimitive("phone").getAsString())
                .name(!contactObject.has("name") ? null : contactObject.getAsJsonPrimitive("name").getAsString())
                .email(!contactObject.has("email") ? null : contactObject.getAsJsonPrimitive("email").getAsString())
                .build();
    }

    @NotNull
    private static List<MolecularMatchTrialsTags> createMolecularMatchTrialsTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTrialsTags> tagsList = Lists.newArrayList();
        for (JsonElement tag : arrayTags) {
            Set<String> keysTags = tag.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found " + keysTags.size() + " ein molecular match trials tags rather than the expected "
                        + EXPECTED_MOLECULARMATCH_TRAILS_TAGS_ELEMENT_SIZES);
                LOGGER.warn(keysTags);
            }
            tagsList.add(ImmutableMolecularMatchTrialsTags.builder()
                    .facet(tag.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .compositeKey(tag.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(tag.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(tag.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tag.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .custom(tag.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .priority(tag.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .alias(!tag.getAsJsonObject().has("alias") ? null : tag.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .manualSuppress(!tag.getAsJsonObject().has("manualSuppress")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedBy(!tag.getAsJsonObject().has("generatedBy")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .generatedByTerm(!tag.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .id(!tag.getAsJsonObject().has("id") ? null : tag.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .manualPriority(!tag.getAsJsonObject().has("manualPriority")
                            ? null
                            : tag.getAsJsonObject().getAsJsonPrimitive("manualPriority").getAsString())
                    .build());
        }
        return tagsList;
    }

    @NotNull
    private static MolecularMatch createMolecularMatch(@NotNull JsonObject objectMolecularMatch) {
        return ImmutableMolecularMatch.builder()
                .criteriaUnmet(createCriteriaUnmet(objectMolecularMatch.getAsJsonArray("criteriaUnmet")))
                .prevalence(createPrevalence(objectMolecularMatch.getAsJsonArray("prevalence")))
                .score(objectMolecularMatch.getAsJsonPrimitive("_score").getAsString())
                .autoGenerateNarrative(objectMolecularMatch.getAsJsonPrimitive("autoGenerateNarrative").getAsString())
                .mutations(createMutations(objectMolecularMatch.getAsJsonArray("mutations")))
                .sources(createSource(objectMolecularMatch.getAsJsonArray("sources")))
                .clinicalSignificance(objectMolecularMatch.getAsJsonPrimitive("clinicalSignificance").getAsString())
                .id(objectMolecularMatch.getAsJsonPrimitive("id").getAsString())
                .includeCondition0(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeCondition0")))
                .includeCondition1(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeCondition1")))
                .uniqueKey(objectMolecularMatch.getAsJsonPrimitive("uniqueKey").getAsString())
                .civicValue(objectMolecularMatch.getAsJsonPrimitive("civic").getAsString())
                .hashKey(objectMolecularMatch.getAsJsonPrimitive("hashKey").getAsString())
                .regulatoryBodyApproved(objectMolecularMatch.getAsJsonPrimitive("regulatoryBodyApproved").getAsString())
                .version(objectMolecularMatch.getAsJsonPrimitive("version").getAsString())
                .includeMutation1(!objectMolecularMatch.has("includeMutation1")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeMutation1")))
                .includeMutation0(!objectMolecularMatch.has("includeMutation0")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeMutation0")))
                .guidelineBody(!objectMolecularMatch.has("guidelineBody")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("guidelineBody").getAsString())
                .regulatoryBody(objectMolecularMatch.getAsJsonPrimitive("regulatoryBody").getAsString())
                .customer(objectMolecularMatch.getAsJsonPrimitive("customer").getAsString())
                .direction(objectMolecularMatch.getAsJsonPrimitive("direction").getAsString())
                .ampcap(objectMolecularMatch.getAsJsonPrimitive("ampcap").getAsString())
                .asts(createAst(objectMolecularMatch.getAsJsonObject("ast")))
                .variantInfo(createVariantInfo(objectMolecularMatch.getAsJsonArray("variantInfo")))
                .guidelineVersion(!objectMolecularMatch.has("guidelineVersion")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("guidelineVersion").getAsString())
                .institution(!objectMolecularMatch.has("institution")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("institution")))
                .tier(objectMolecularMatch.getAsJsonPrimitive("tier").getAsString())
                .tierExplanation(createTierExplanation(objectMolecularMatch.getAsJsonArray("tierExplanation")))
                .mvld(objectMolecularMatch.getAsJsonPrimitive("mvld").getAsString())
                .tags(createTags(objectMolecularMatch.getAsJsonArray("tags")))
                .criteriaMet(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("criteriaMet")))
                .biomarkerClass(objectMolecularMatch.getAsJsonPrimitive("biomarkerClass").getAsString())
                .classification(createClassification(objectMolecularMatch.getAsJsonArray("classifications")))
                .includeDrug1(!objectMolecularMatch.has("includeDrug1")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeDrug1")))
                .includeStage0(!objectMolecularMatch.has("includeStage0")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeStage0")))
                .therapeuticContext(createTherapeuticContext(objectMolecularMatch.getAsJsonArray("therapeuticContext")))
                .sixtier(objectMolecularMatch.getAsJsonPrimitive("sixtier").getAsString())
                .noTherapyAvailable(!objectMolecularMatch.has("noTherapyAvailable")
                        ? null
                        : objectMolecularMatch.getAsJsonPrimitive("noTherapyAvailable").getAsString())
                .external_id(!objectMolecularMatch.has("external_id")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("external_id")))
                .narrative(objectMolecularMatch.getAsJsonPrimitive("narrative").getAsString())
                .expression(objectMolecularMatch.getAsJsonPrimitive("expression").getAsString())
                .includeGene0(!objectMolecularMatch.has("includeDrug0")
                        ? null
                        : jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeGene0")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchTherapeuticContext> createTherapeuticContext(@NotNull JsonArray arrayTherapeuticContext) {
        List<MolecularMatchTherapeuticContext> therapeuticContextList = Lists.newArrayList();
        for (JsonElement therapeuticContext : arrayTherapeuticContext) {
            Set<String> keysTherapeuticContext = therapeuticContext.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES.contains(keysTherapeuticContext.size())) {
                LOGGER.warn("Found " + keysTherapeuticContext.size() + " in molecular match therapeutic context rather than the expected "
                        + EXPECTED_MOLECULARMATCH_THERAPEUTIC_CONTEXT_SIZES);
                LOGGER.warn(keysTherapeuticContext);
            }

            therapeuticContextList.add(ImmutableMolecularMatchTherapeuticContext.builder()
                    .facet(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .suppress(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .valid(!therapeuticContext.getAsJsonObject().has("valid")
                            ? null
                            : therapeuticContext.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .name(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return therapeuticContextList;
    }

    @NotNull
    private static List<MolecularMatchClassification> createClassification(@NotNull JsonArray objectClassifications) {

        List<MolecularMatchClassification> classificationList = Lists.newArrayList();
        for (JsonElement classification : objectClassifications) {
            Set<String> keysClassification = classification.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES.contains(keysClassification.size())) {
                LOGGER.warn("Found " + keysClassification.size() + " in molecular match classification rather than the expected "
                        + EXPECTED_MOLECULARMATCH_CLASSIFICATION_SIZES);
                LOGGER.warn(keysClassification);
            }

            classificationList.add(ImmutableMolecularMatchClassification.builder()
                    .end(!classification.getAsJsonObject().has("End")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("End")))
                    .classification(classification.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .classificationOverride(
                            !classification.getAsJsonObject().has("classificationOverride") || classification.getAsJsonObject()
                                    .get("classificationOverride")
                                    .isJsonNull()
                                    ? null
                                    : classification.getAsJsonObject().getAsJsonPrimitive("classificationOverride").getAsString())
                    .start(!classification.getAsJsonObject().has("Start")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Start")))
                    .chr(!classification.getAsJsonObject().has("Chr")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Chr")))
                    .geneSymbol(!classification.getAsJsonObject().has("geneSymbol")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(!classification.getAsJsonObject().has("pathology")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("pathology")))
                    .ref(!classification.getAsJsonObject().has("Ref")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Ref")))
                    .description(!classification.getAsJsonObject().has("description")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .priority(!classification.getAsJsonObject().has("priority")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .NucleotideChange(!classification.getAsJsonObject().has("NucleotideChange")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("NucleotideChange")))
                    .parents(!classification.getAsJsonObject().has("parents")
                            ? null
                            : createParents(classification.getAsJsonObject().getAsJsonArray("parents")))
                    .expandGeneSearch(!classification.getAsJsonObject().has("expandGeneSearch")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("expandGeneSearch").getAsString())
                    .drugsExperimentalCount(!classification.getAsJsonObject().has("drugsExperimentalCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsExperimentalCount").getAsString())
                    .exon(!classification.getAsJsonObject().has("Exon")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Exon")))
                    .drugsApprovedOffLabelCount(!classification.getAsJsonObject().has("drugsApprovedOffLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOffLabelCount").getAsString())
                    .exonicFunc(!classification.getAsJsonObject().has("ExonicFunc")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("ExonicFunc")))
                    .popFreqMax(!classification.getAsJsonObject().has("PopFreqMax")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("PopFreqMax")))
                    .copyNumberType(!classification.getAsJsonObject().has("copyNumberType") || classification.getAsJsonObject()
                            .get("copyNumberType")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("copyNumberType").getAsString())
                    .publicationCount(!classification.getAsJsonObject().has("publicationCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("publicationCount").getAsString())
                    .transcript(!classification.getAsJsonObject().has("transcript") || classification.getAsJsonObject()
                            .get("transcript")
                            .isJsonNull() ? null : classification.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .dbSNP(!classification.getAsJsonObject().has("dbSNP")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("dbSNP")))
                    .alt(!classification.getAsJsonObject().has("Alt")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Alt")))
                    .name(!classification.getAsJsonObject().has("name")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .rootTerm(!classification.getAsJsonObject().has("rootTerm")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("rootTerm").getAsString())
                    .sources(!classification.getAsJsonObject().has("sources")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("sources")))
                    .drugsApprovedOnLabelCount(!classification.getAsJsonObject().has("drugsApprovedOnLabelCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOnLabelCount").getAsString())
                    .trialCount(!classification.getAsJsonObject().has("trialCount")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("trialCount").getAsString())
                    .alias(!classification.getAsJsonObject().has("alias")
                            ? null
                            : classification.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .COSMIC_ID(!classification.getAsJsonObject().has("COSMIC_ID")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("COSMIC_ID")))
                    .transcripts(!classification.getAsJsonObject().has("transcripts")
                            ? null
                            : jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("transcripts")))
                    .build());

        }
        return classificationList;

    }

    @NotNull
    private static List<MolecularMatchParents> createParents(@NotNull JsonArray arrayParents) {
        List<MolecularMatchParents> parentsList = Lists.newArrayList();
        for (JsonElement parents : arrayParents) {
            Set<String> keysParents = parents.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PARENTS_SIZES.contains(keysParents.size())) {
                LOGGER.warn("Found " + keysParents.size() + " in molecular match parents rather than the expected "
                        + EXPECTED_MOLECULARMATCH_PARENTS_SIZES);
                LOGGER.warn(keysParents);
            }
            parentsList.add(ImmutableMolecularMatchParents.builder()
                    .transcripts(jsonArrayToStringList(parents.getAsJsonObject().getAsJsonArray("transcripts")))
                    .type(!parents.getAsJsonObject().has("type") || parents.getAsJsonObject().get("type").isJsonNull()
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .name(parents.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .actionableParent(!parents.getAsJsonObject().has("actionableParent")
                            ? null
                            : parents.getAsJsonObject().getAsJsonPrimitive("actionableParent").getAsString())
                    .build());
        }

        return parentsList;
    }

    @NotNull
    private static List<MolecularMatchTags> createTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTags> tagsList = Lists.newArrayList();
        for (JsonElement tags : arrayTags) {
            Set<String> keysTags = tags.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TAGS_SIZES.contains(keysTags.size())) {
                LOGGER.warn("Found " + keysTags.size() + " in molecular match tags rather than the expected "
                        + EXPECTED_MOLECULARMATCH_TAGS_SIZES);
                LOGGER.warn(keysTags);
            }

            tagsList.add(ImmutableMolecularMatchTags.builder()
                    .priority(tags.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(!tags.getAsJsonObject().has("compositeKey")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(!tags.getAsJsonObject().has("suppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(!tags.getAsJsonObject().has("filterType")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tags.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!tags.getAsJsonObject().has("primary")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(tags.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!tags.getAsJsonObject().has("valid") ? null : tags.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!tags.getAsJsonObject().has("custom")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .isNew(!tags.getAsJsonObject().has("isNew") ? null : tags.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!tags.getAsJsonObject().has("generatedBy")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!tags.getAsJsonObject().has("manualSuppress")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!tags.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .transcript(!tags.getAsJsonObject().has("transcript")
                            ? null
                            : tags.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return tagsList;
    }

    @NotNull
    private static List<MolecularMatchTierExplanation> createTierExplanation(@NotNull JsonArray arrarTierExplanation) {
        List<MolecularMatchTierExplanation> tierExplanationList = Lists.newArrayList();
        for (JsonElement tierExplanation : arrarTierExplanation) {
            Set<String> keysTierExplanation = tierExplanation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES.contains(keysTierExplanation.size())) {
                LOGGER.warn("Found " + keysTierExplanation.size() + " in molecular match tier explanation rather than the expected "
                        + EXPECTED_MOLECULARMATCH_TIER_EXPLANATION_SIZES);
                LOGGER.warn(keysTierExplanation);
            }

            tierExplanationList.add(ImmutableMolecularMatchTierExplanation.builder()
                    .tier(tierExplanation.getAsJsonObject().getAsJsonPrimitive("tier").getAsString())
                    .step(tierExplanation.getAsJsonObject().getAsJsonPrimitive("step").getAsString())
                    .message(tierExplanation.getAsJsonObject().getAsJsonPrimitive("message").getAsString())
                    .success(tierExplanation.getAsJsonObject().getAsJsonPrimitive("success").getAsString())
                    .build());
        }
        return tierExplanationList;
    }

    @NotNull
    private static List<MolecularMatchVariantInfo> createVariantInfo(@NotNull JsonArray arrayVariantInfo) {
        List<MolecularMatchVariantInfo> variantInfoList = Lists.newArrayList();

        for (JsonElement variantInfo : arrayVariantInfo) {
            Set<String> keysVariantInfo = variantInfo.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES.contains(keysVariantInfo.size())) {
                LOGGER.warn("Found " + keysVariantInfo.size() + " in molecular match variant info rather than the expected "
                        + EXPECTED_MOLECULARMATCH_VARIANTINFO_SIZES);
                LOGGER.warn(keysVariantInfo);
            }
            variantInfoList.add(ImmutableMolecularMatchVariantInfo.builder()
                    .classification(variantInfo.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .name(variantInfo.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .consequences(jsonArrayToStringList(variantInfo.getAsJsonObject().getAsJsonArray("consequences")))
                    .fusions(createFusions(variantInfo.getAsJsonObject().getAsJsonArray("fusions")))
                    .locations(createLocations(variantInfo.getAsJsonObject().getAsJsonArray("locations")))
                    .geneFusionPartner(variantInfo.getAsJsonObject().getAsJsonPrimitive("geneFusionPartner").getAsString())
                    .COSMIC_ID(variantInfo.getAsJsonObject().get("COSMIC_ID").isJsonNull()
                            ? null
                            : variantInfo.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .gene(variantInfo.getAsJsonObject().getAsJsonPrimitive("gene").getAsString())
                    .transcript(variantInfo.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .popFreqMax(variantInfo.getAsJsonObject().getAsJsonPrimitive("popFreqMax").getAsString())
                    .build());
        }
        return variantInfoList;
    }

    @NotNull
    private static List<MolecularMatchFusions> createFusions(@NotNull JsonArray arrayFusions) {
        List<MolecularMatchFusions> fusionsList = Lists.newArrayList();

        for (JsonElement fusions : arrayFusions) {
            Set<String> keysFusions = fusions.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONS_SIZES.contains(keysFusions.size())) {
                LOGGER.warn("Found " + keysFusions.size() + " in molecular match variant info rather than the expected "
                        + EXPECTED_MOLECULARMATCH_FUSIONS_SIZES);
                LOGGER.warn(keysFusions);
            }
            fusionsList.add(ImmutableMolecularMatchFusions.builder()
                    .referenceGenome(fusions.getAsJsonObject().get("referenceGenome").getAsString())
                    .LBPWREP(fusions.getAsJsonObject().get("LBPWREP").getAsString())
                    .RBPWREP(fusions.getAsJsonObject().get("RBPWREP").getAsString())
                    .exonNumber(fusions.getAsJsonObject().get("exonNumber").getAsString())
                    .chr(fusions.getAsJsonObject().get("chr").getAsString())
                    .RBPWLEP(fusions.getAsJsonObject().get("RBPWLEP").getAsString())
                    .intronNumber(fusions.getAsJsonObject().get("intronNumber").getAsString())
                    .LBPWLEP(fusions.getAsJsonObject().get("LBPWLEP").getAsString())
                    .build());
        }
        return fusionsList;
    }

    @NotNull
    private static List<MolecularMatchLocations> createLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchLocations> locationsList = Lists.newArrayList();
        for (JsonElement locations : arrayLocations) {
            Set<String> keysLocations = locations.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES.contains(keysLocations.size())) {
                LOGGER.warn("Found " + keysLocations.size() + " in molecular match locations rather than the expected "
                        + EXPECTED_MOLECULARMATCH_LOCATIONS_SIZES);
                LOGGER.warn(keysLocations);
            }

            locationsList.add(ImmutableMolecularMatchLocations.builder()
                    .aminoAcidChange(!locations.getAsJsonObject().has("amino_acid_change")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .intronNumber(!locations.getAsJsonObject().has("intronNumber")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(!locations.getAsJsonObject().has("exonNumber") ? null : createArrayExonNumber(locations))
                    .stop(locations.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(locations.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(locations.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(!locations.getAsJsonObject().has("strand")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .alt(!locations.getAsJsonObject().has("alt")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .referenceGenome(!locations.getAsJsonObject().has("referenceGenome")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(!locations.getAsJsonObject().has("ref")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .cdna(!locations.getAsJsonObject().has("cdna")
                            ? null
                            : locations.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return locationsList;
    }

    @NotNull
    private static List<String> createArrayExonNumber(@NotNull JsonElement elementExonNumber) {
        if (elementExonNumber.getAsJsonObject().get("exonNumber").isJsonArray()) {
            return jsonArrayToStringList(elementExonNumber.getAsJsonObject().getAsJsonArray("exonNumber"));
        } else {
            return Arrays.asList(elementExonNumber.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString());
        }
    }

    @NotNull
    private static MolecularMatchAst createAst(@NotNull JsonObject objectAst) {
        Set<String> keysAst = objectAst.keySet();
        if (!EXPECTED_MOLECULARMATCH_AST_SIZES.contains(keysAst.size())) {
            LOGGER.warn(
                    "Found " + keysAst.size() + " in molecular match ast rather than the expected " + EXPECTED_MOLECULARMATCH_AST_SIZES);
            LOGGER.warn(keysAst);
        }
        return ImmutableMolecularMatchAst.builder()
                .raw(objectAst.get("raw") == null ? null : objectAst.getAsJsonPrimitive("raw").getAsString())
                .value(!objectAst.has("value") ? null : objectAst.getAsJsonPrimitive("value").getAsString())
                .operator(!objectAst.has("operator") ? null : objectAst.getAsJsonPrimitive("operator").getAsString())
                .right(!objectAst.has("right") ? null : createRight(objectAst.getAsJsonObject("right")))
                .type(objectAst.getAsJsonPrimitive("type").getAsString())
                .left(!objectAst.has("left") ? null : createLeft(objectAst.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstLeft createLeft(@NotNull JsonObject objectLeft) {
        Set<String> keysLeft = objectLeft.keySet();
        if (!EXPECTED_MOLECULARMATCH_LEFT_SIZES.contains(keysLeft.size())) {
            LOGGER.warn("Found " + keysLeft.size() + " in molecular match ast left rather than the expected "
                    + EXPECTED_MOLECULARMATCH_LEFT_SIZES);
            LOGGER.warn(keysLeft);
        }
        return ImmutableMolecularMatchAstLeft.builder()
                .operator(!objectLeft.has("operator") ? null : objectLeft.getAsJsonPrimitive("operator").getAsString())
                .raw(!objectLeft.has("raw") ? null : objectLeft.getAsJsonPrimitive("raw").getAsString())
                .type(objectLeft.getAsJsonPrimitive("type").getAsString())
                .value(!objectLeft.has("value") ? null : objectLeft.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRight createRight(@NotNull JsonObject objectRight) {
        Set<String> keysRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_SIZES.contains(keysRight.size())) {
            LOGGER.warn("Found " + keysRight.size() + " in molecular match ast right rather than the expected "
                    + EXPECTED_MOLECULARMATCH_RIGHT_SIZES);
            LOGGER.warn(keysRight);
        }
        return ImmutableMolecularMatchAstRight.builder()
                .operator(!objectRight.has("operator") ? null : objectRight.getAsJsonPrimitive("operator").getAsString())
                .left(!objectRight.has("left") ? null : createRightLeft(objectRight.getAsJsonObject("left")))
                .right(!objectRight.has("right") ? null : createRightRight(objectRight.getAsJsonObject("right")))
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();

    }

    @NotNull
    private static MolecularMatchAstRightRight createRightRight(@NotNull JsonObject objectRight) {
        Set<String> keysRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES.contains(keysRight.size())) {
            LOGGER.warn("Found " + keysRight.size() + " in molecular match right right rather than the expected "
                    + EXPECTED_MOLECULARMATCH_RIGHT_RIGHT_SIZES);
            LOGGER.warn(keysRight);
        }
        return ImmutableMolecularMatchAstRightRight.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRightLeft createRightLeft(@NotNull JsonObject objectRight) {
        Set<String> keysRight = objectRight.keySet();
        if (!EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES.contains(keysRight.size())) {
            LOGGER.warn("Found " + keysRight.size() + " in molecular match right left rather than the expected "
                    + EXPECTED_MOLECULARMATCH_RIGHT_LEFT_SIZES);
            LOGGER.warn(keysRight);
        }
        return ImmutableMolecularMatchAstRightLeft.builder()
                .raw(!objectRight.has("raw") ? null : objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(!objectRight.has("value") ? null : objectRight.getAsJsonPrimitive("value").getAsString())
                .left(!objectRight.has("left") ? null : createRight(objectRight.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchSource> createSource(@NotNull JsonArray arraySources) {
        List<MolecularMatchSource> sourcesList = Lists.newArrayList();
        for (JsonElement source : arraySources) {
            Set<String> keysSource = source.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_SOURCE_SIZES.contains(keysSource.size())) {
                LOGGER.warn("Found " + keysSource.size() + " in molecular match source rather than the expected "
                        + EXPECTED_MOLECULARMATCH_SOURCE_SIZES);
                LOGGER.warn(keysSource);
            }

            sourcesList.add(ImmutableMolecularMatchSource.builder()
                    .name(source.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .suppress(source.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .pubId(source.getAsJsonObject().getAsJsonPrimitive("pubId").getAsString())
                    .subType(!source.getAsJsonObject().has("subType")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("subType").getAsString())
                    .valid(source.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .link(source.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .year(source.getAsJsonObject().getAsJsonPrimitive("year").getAsString())
                    .trialId(!source.getAsJsonObject().has("trialId")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialId").getAsString())
                    .type(source.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(source.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .institution(!source.getAsJsonObject().has("institution")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("institution").getAsString())
                    .trialPhase(!source.getAsJsonObject().has("trialPhase")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trialPhase").getAsString())
                    .functionalConsequence(!source.getAsJsonObject().has("functionalConsequence")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("functionalConsequence").getAsString())
                    .trustRating(!source.getAsJsonObject().has("trustRating")
                            ? null
                            : source.getAsJsonObject().getAsJsonPrimitive("trustRating").getAsString())
                    .build());
        }
        return sourcesList;
    }

    @NotNull
    private static List<MolecularMatchCriteriaUnmet> createCriteriaUnmet(@NotNull JsonArray arrayCriteriaUnmet) {
        List<MolecularMatchCriteriaUnmet> criteriaUnmetList = Lists.newArrayList();

        for (JsonElement criteriaUnmet : arrayCriteriaUnmet) {
            Set<String> keysCriteriaUnmet = criteriaUnmet.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES.contains(keysCriteriaUnmet.size())) {
                LOGGER.warn("Found " + keysCriteriaUnmet.size() + " in molecular match criteria unmet rather than the expected "
                        + EXPECTED_MOLECULARMATCH_CRITERIA_UNMET_SIZES);
                LOGGER.warn(keysCriteriaUnmet);
            }

            criteriaUnmetList.add(ImmutableMolecularMatchCriteriaUnmet.builder()
                    .priority(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .isNew(!criteriaUnmet.getAsJsonObject().has("isNew")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("isNew").getAsString())
                    .generatedBy(!criteriaUnmet.getAsJsonObject().has("generatedBy")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedBy").getAsString())
                    .manualSuppress(!criteriaUnmet.getAsJsonObject().has("manualSuppress")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("manualSuppress").getAsString())
                    .generatedByTerm(!criteriaUnmet.getAsJsonObject().has("generatedByTerm")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("generatedByTerm").getAsString())
                    .suppress(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(!criteriaUnmet.getAsJsonObject().has("primary")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(!criteriaUnmet.getAsJsonObject().has("valid")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(!criteriaUnmet.getAsJsonObject().has("custom")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .transcript(!criteriaUnmet.getAsJsonObject().has("transcript")
                            ? null
                            : criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .build());
        }
        return criteriaUnmetList;
    }

    @NotNull
    private static List<MolecularMatchPrevalence> createPrevalence(@NotNull JsonArray arrayPrevelance) {
        List<MolecularMatchPrevalence> prevalenceList = Lists.newArrayList();

        for (JsonElement prevelance : arrayPrevelance) {
            Set<String> keysPrevalence = prevelance.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_PREVELANCE_SIZES.contains(keysPrevalence.size())) {
                LOGGER.warn("Found " + keysPrevalence.size() + " in molecular match prevalence rather than the expected "
                        + EXPECTED_MOLECULARMATCH_PREVELANCE_SIZES);
                LOGGER.warn(keysPrevalence);
            }

            prevalenceList.add(ImmutableMolecularMatchPrevalence.builder()
                    .count(prevelance.getAsJsonObject().getAsJsonPrimitive("count").getAsString())
                    .percent(prevelance.getAsJsonObject().getAsJsonPrimitive("percent").getAsString())
                    .studyId(prevelance.getAsJsonObject().getAsJsonPrimitive("studyId").getAsString())
                    .samples(prevelance.getAsJsonObject().getAsJsonPrimitive("samples").getAsString())
                    .molecular(!prevelance.getAsJsonObject().has("molecular")
                            ? null
                            : prevelance.getAsJsonObject().getAsJsonPrimitive("molecular").getAsString())
                    .condition(!prevelance.getAsJsonObject().has("condition")
                            ? null
                            : prevelance.getAsJsonObject().getAsJsonPrimitive("condition").getAsString())
                    .build());
        }
        return prevalenceList;
    }

    @NotNull
    private static List<MolecularMatchMutations> createMutations(@NotNull JsonArray arrayMutations) {
        List<MolecularMatchMutations> mutationList = Lists.newArrayList();

        for (JsonElement mutation : arrayMutations) {
            Set<String> keysMutations = mutation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES.contains(keysMutations.size())) {
                LOGGER.warn("Found " + keysMutations.size() + " in molecular match mutations rather than the expected "
                        + EXPECTED_MOLECULARMATCH_MUTATIONS_SIZES);
                LOGGER.warn(keysMutations);
            }

            mutationList.add(ImmutableMolecularMatchMutations.builder()
                    .transcriptConsequence(mutation.getAsJsonObject().get("transcriptConsequence") == null
                            ? null
                            : createTranscriptConsequence(mutation.getAsJsonObject().getAsJsonArray("transcriptConsequence")))
                    .longestTranscript(!mutation.getAsJsonObject().has("longestTranscript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("longestTranscript").getAsString())
                    .parents(createParents(mutation.getAsJsonObject().getAsJsonArray("parents")))
                    .wgsaData(!mutation.getAsJsonObject().has("wgsaData")
                            ? null
                            : createWgsaData(mutation.getAsJsonObject().getAsJsonObject("wgsaData")))
                    .wgsaMap(!mutation.getAsJsonObject().has("wgsaMap")
                            ? null
                            : createWgsaMap(mutation.getAsJsonObject().getAsJsonArray("wgsaMap")))
                    .exonsInfo(!mutation.getAsJsonObject().has("exonsInfo")
                            ? null
                            : createExonsInfo(mutation.getAsJsonObject().getAsJsonObject("exonsInfo")))
                    .fusionData(!mutation.getAsJsonObject().has("fusionData")
                            ? null
                            : createFusionData(mutation.getAsJsonObject().getAsJsonArray("fusionData")))
                    .transcriptRecognized(!mutation.getAsJsonObject().has("transcriptRecognized")
                            ? null
                            : mutation.getAsJsonObject().get("transcriptRecognized").getAsString())
                    .description(mutation.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .mutationType(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("mutation_type")))
                    .src(mutation.getAsJsonObject().getAsJsonPrimitive("_src").getAsString())
                    .sources(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("sources")))
                    .synonyms(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("synonyms")))
                    .gRch37Location(createGRCH37Location(mutation.getAsJsonObject().getAsJsonArray("GRCh37_location")))
                    .uniprotTranscript(!mutation.getAsJsonObject().has("uniprotTranscript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("uniprotTranscript").getAsString())
                    .geneSymbol(mutation.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("pathology")))
                    .transcript(!mutation.getAsJsonObject().has("transcript")
                            ? null
                            : mutation.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .id(mutation.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .cDNA(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("cdna")))
                    .name(mutation.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());

        }
        return mutationList;
    }

    @NotNull
    private static List<MolecularMatchFusionData> createFusionData(@NotNull JsonArray arrayFusionData) {
        List<MolecularMatchFusionData> fusionDataList = Lists.newArrayList();
        for (JsonElement fusiondata : arrayFusionData.getAsJsonArray()) {
            Set<String> keysFusionData = fusiondata.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES.contains(keysFusionData.size())) {
                LOGGER.warn("Found " + keysFusionData.size() + " in molecular match fusion data rather than the expected "
                        + EXPECTED_MOLECULARMATCH_FUSIONDATA_SIZES);
                LOGGER.warn(keysFusionData);
            }

            fusionDataList.add(ImmutableMolecularMatchFusionData.builder()
                    .Bgreg(!fusiondata.getAsJsonObject().has("Bgreg")
                            ? null
                            : createBreg(fusiondata.getAsJsonObject().getAsJsonArray("Bgreg")))
                    .Bchr(!fusiondata.getAsJsonObject().has("Bchr")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Bchr")))
                    .synonym(!fusiondata.getAsJsonObject().has("synonym")
                            ? null
                            : fusiondata.getAsJsonObject().getAsJsonPrimitive("synonym").getAsString())
                    .Agene(!fusiondata.getAsJsonObject().has("Agene")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Agene")))
                    .Btx(!fusiondata.getAsJsonObject().has("Btx")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Btx")))
                    .Achr(!fusiondata.getAsJsonObject().has("Achr")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Achr")))
                    .ins(!fusiondata.getAsJsonObject().has("ins")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("ins")))
                    .source(!fusiondata.getAsJsonObject().has("source")
                            ? null
                            : fusiondata.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .Agreg(!fusiondata.getAsJsonObject().has("Agreg")
                            ? null
                            : createAreg(fusiondata.getAsJsonObject().getAsJsonArray("Agreg")))
                    .Bgene(!fusiondata.getAsJsonObject().has("Bgene")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Bgene")))
                    .Acoord(!fusiondata.getAsJsonObject().has("Acoord")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Acoord")))
                    .Bori(!fusiondata.getAsJsonObject().has("Bori")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Bori")))
                    .Aband(!fusiondata.getAsJsonObject().has("Aband")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Aband")))
                    .Bband(!fusiondata.getAsJsonObject().has("Bband")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Bband")))
                    .Aori(!fusiondata.getAsJsonObject().has("Aori")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Aori")))
                    .Atx(!fusiondata.getAsJsonObject().has("Atx")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Atx")))
                    .Bcoord(!fusiondata.getAsJsonObject().has("Bcoord")
                            ? null
                            : jsonArrayToStringList(fusiondata.getAsJsonObject().getAsJsonArray("Bcoord")))
                    .Paper(!fusiondata.getAsJsonObject().has("Paper")
                            ? null
                            : fusiondata.getAsJsonObject().getAsJsonPrimitive("Paper").getAsString())
                    .build());
        }
        return fusionDataList;
    }

    @NotNull
    private static List<MolecularMatchAreg> createAreg(@NotNull JsonArray arrayAreg) {
        List<MolecularMatchAreg> AregList = Lists.newArrayList();
        for (JsonElement areg : arrayAreg.getAsJsonArray()) {
            Set<String> keysAreg = areg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_AREG_SIZES.contains(keysAreg.size())) {
                LOGGER.warn("Found " + keysAreg.size() + " in molecular match areg rather than the expected "
                        + EXPECTED_MOLECULARMATCH_AREG_SIZES);
                LOGGER.warn(keysAreg);
            }
            AregList.add(ImmutableMolecularMatchAreg.builder()
                    .type(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(areg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .build());
        }
        return AregList;
    }

    @NotNull
    private static List<MolecularMatchBreg> createBreg(@NotNull JsonArray arrayBreg) {
        List<MolecularMatchBreg> bregList = Lists.newArrayList();
        for (JsonElement breg : arrayBreg.getAsJsonArray()) {
            Set<String> keysBreg = breg.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_BREG_SIZES.contains(keysBreg.size())) {
                LOGGER.warn("Found " + keysBreg.size() + " in molecular match breg rather than the expected "
                        + EXPECTED_MOLECULARMATCH_BREG_SIZES);
                LOGGER.warn(keysBreg);
            }
            bregList.add(ImmutableMolecularMatchBreg.builder()
                    .type(breg.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .num(breg.getAsJsonObject().getAsJsonPrimitive("num").getAsString())
                    .build());
        }
        return bregList;
    }

    @NotNull
    private static MolecularMatchExonInfo createExonsInfo(@NotNull JsonObject objectExonsInfo) {
        Set<String> keysExonsInfo = objectExonsInfo.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES.contains(keysExonsInfo.size())) {
            LOGGER.warn("Found " + keysExonsInfo.size() + " in molecular match exons info rather than the expected "
                    + EXPECTED_MOLECULARMATCH_EXONSINFO_SIZES);
            LOGGER.warn(keysExonsInfo);
        }
        return ImmutableMolecularMatchExonInfo.builder()
                .exonBoundaries(createExonBoundries(objectExonsInfo.getAsJsonObject("exonBoundaries")))
                .txStart(!objectExonsInfo.has("txStart") ? null : objectExonsInfo.getAsJsonPrimitive("txStart").getAsString())
                .cdsEnd(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .chr(objectExonsInfo.getAsJsonPrimitive("chr").getAsString())
                .cdsStart(!objectExonsInfo.has("cdsEnd") ? null : objectExonsInfo.getAsJsonPrimitive("cdsEnd").getAsString())
                .transcript(objectExonsInfo.getAsJsonPrimitive("transcript").getAsString())
                .txEnd(!objectExonsInfo.has("txEnd") ? null : objectExonsInfo.getAsJsonPrimitive("txEnd").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchPositions createMolecularPositions(@NotNull JsonObject objectPositions) {
        Set<String> keysPositions = objectPositions.keySet();
        if (!EXPECTED_MOLECULARMATCH_POSITIONS_SIZES.contains(keysPositions.size())) {
            LOGGER.warn("Found " + keysPositions.size() + " in molecular match positions rather than the expected "
                    + EXPECTED_MOLECULARMATCH_POSITIONS_SIZES);
            LOGGER.warn(keysPositions);
        }
        return ImmutableMolecularMatchPositions.builder()
                .start(objectPositions.getAsJsonPrimitive("start").getAsString())
                .stop(objectPositions.getAsJsonPrimitive("stop").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchExonsBoundries createExonBoundries(@NotNull JsonObject objectExonsBoundries) {
        Set<String> keysExonsBoundries = objectExonsBoundries.keySet();
        if (!EXPECTED_MOLECULARMATCH_EXONSBOUNDRIES_SIZES.contains(keysExonsBoundries.size())) {
            LOGGER.warn("Found " + keysExonsBoundries.size() + " in molecular match boundries rather than the expected "
                    + EXPECTED_MOLECULARMATCH_EXONSBOUNDRIES_SIZES);
            LOGGER.warn(keysExonsBoundries);
        }

        return ImmutableMolecularMatchExonsBoundries.builder()
                .exon1(!objectExonsBoundries.has("1") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("1")))
                .exon2(!objectExonsBoundries.has("2") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("2")))
                .exon3(!objectExonsBoundries.has("3") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("3")))
                .exon4(!objectExonsBoundries.has("4") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("4")))
                .exon5(!objectExonsBoundries.has("5") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("5")))
                .exon6(!objectExonsBoundries.has("6") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("6")))
                .exon7(!objectExonsBoundries.has("7") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("7")))
                .exon8(!objectExonsBoundries.has("8") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("8")))
                .exon9(!objectExonsBoundries.has("9") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("9")))
                .exon10(!objectExonsBoundries.has("10") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("10")))
                .exon11(!objectExonsBoundries.has("11") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("11")))
                .exon12(!objectExonsBoundries.has("12") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("12")))
                .exon13(!objectExonsBoundries.has("13") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("13")))
                .exon14(!objectExonsBoundries.has("14") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("14")))
                .exon15(!objectExonsBoundries.has("15") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("15")))
                .exon16(!objectExonsBoundries.has("16") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("16")))
                .exon17(!objectExonsBoundries.has("17") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("17")))
                .exon18(!objectExonsBoundries.has("18") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("18")))
                .exon19(!objectExonsBoundries.has("19") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("19")))
                .exon20(!objectExonsBoundries.has("20") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("20")))
                .exon21(!objectExonsBoundries.has("21") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("21")))
                .exon22(!objectExonsBoundries.has("22") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("22")))
                .exon23(!objectExonsBoundries.has("23") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("23")))
                .exon24(!objectExonsBoundries.has("24") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("24")))
                .exon25(!objectExonsBoundries.has("25") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("25")))
                .exon26(!objectExonsBoundries.has("26") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("26")))
                .exon27(!objectExonsBoundries.has("27") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("27")))
                .exon28(!objectExonsBoundries.has("28") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("28")))
                .exon29(!objectExonsBoundries.has("29") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("29")))
                .exon30(!objectExonsBoundries.has("30") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("30")))
                .exon31(!objectExonsBoundries.has("31") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("31")))
                .exon32(!objectExonsBoundries.has("32") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("32")))
                .exon33(!objectExonsBoundries.has("33") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("33")))
                .exon34(!objectExonsBoundries.has("34") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("34")))
                .exon35(!objectExonsBoundries.has("35") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("35")))
                .exon36(!objectExonsBoundries.has("36") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("36")))
                .exon37(!objectExonsBoundries.has("37") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("37")))
                .exon38(!objectExonsBoundries.has("38") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("38")))
                .exon39(!objectExonsBoundries.has("39") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("39")))
                .exon40(!objectExonsBoundries.has("40") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("40")))
                .exon41(!objectExonsBoundries.has("41") ? null : createMolecularPositions(objectExonsBoundries.getAsJsonObject("41")))
                .build();
    }

    @NotNull
    private static List<MolecularMatchWGSaMap> createWgsaMap(@NotNull JsonArray objectWgsaMap) {
        List<MolecularMatchWGSaMap> molecularMatchWGSaMapList = Lists.newArrayList();
        for (JsonElement wgsDataMap : objectWgsaMap.getAsJsonArray()) {
            Set<String> keysWgsaMap = wgsDataMap.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES.contains(keysWgsaMap.size())) {
                LOGGER.warn("Found " + keysWgsaMap.size() + " in molecular match wgsa map rather than the expected "
                        + EXPECTED_MOLECULARMATCH_WGSAMAP_SIZES);
                LOGGER.warn(keysWgsaMap);
            }
            molecularMatchWGSaMapList.add(ImmutableMolecularMatchWGSaMap.builder()
                    .AA(!wgsDataMap.getAsJsonObject().has("AA")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("AA").getAsString())
                    .name(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .GRCh37_Chr_Start_Ref_Alt(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("GRCh37_Chr_Start_Ref_Alt").getAsString())
                    .Synonyms(jsonArrayToStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("Synonyms")))
                    .ProtCoords(jsonArrayToStringList(wgsDataMap.getAsJsonObject().getAsJsonArray("ProtCoords")))
                    .NucleotideChange(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("NucleotideChange").getAsString())
                    .Exon(!wgsDataMap.getAsJsonObject().has("Exon")
                            ? null
                            : wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Exon").getAsString())
                    .Gene(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Gene").getAsString())
                    .Transcript(wgsDataMap.getAsJsonObject().getAsJsonPrimitive("Transcript").getAsString())
                    .build());
        }

        return molecularMatchWGSaMapList;
    }

    @NotNull
    private static List<MolecularMatchWGSadataLocation> createWgsaData(@NotNull JsonObject objectWgsaData) {
        List<MolecularMatchWGSadataLocation> molecularMatchWGSadataLocationList = Lists.newArrayList();
        Set<String> keysWgsaData = objectWgsaData.keySet();
        if (!EXPECTED_MOLECULARMATCH_WGSADATA_SIZES.contains(keysWgsaData.size())) {
            LOGGER.warn("Found " + keysWgsaData.size() + " in molecular match wgsa data rather than the expected "
                    + EXPECTED_MOLECULARMATCH_WGSADATA_SIZES);
            LOGGER.warn(keysWgsaData);
        }

        for (JsonElement wgsDataLocation : objectWgsaData.get("locations").getAsJsonArray()) {
            Set<String> keysWgsaDataLocation = wgsDataLocation.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES.contains(keysWgsaDataLocation.size())) {
                LOGGER.warn("Found " + keysWgsaDataLocation.size() + " in molecular match wgsa data location rather than the expected "
                        + EXPECTED_MOLECULARMATCH_WGSADATA_LOCATION_SIZES);
                LOGGER.warn(keysWgsaDataLocation);
            }

            molecularMatchWGSadataLocationList.add(ImmutableMolecularMatchWGSadataLocation.builder()
                    .ExonicFunc(!wgsDataLocation.getAsJsonObject().has("ExonicFunc")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExonicFunc").getAsString())
                    .dbSNP(!wgsDataLocation.getAsJsonObject().has("dbSNP")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("dbSNP").getAsString())
                    .ClinVar_DIS(!wgsDataLocation.getAsJsonObject().has("ClinVar_DIS")
                            ? null
                            : jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_DIS")))
                    .ClinVar_SIG(!wgsDataLocation.getAsJsonObject().has("ClinVar_SIG")
                            ? null
                            : jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_SIG")))
                    .ClinVar_STATUS(!wgsDataLocation.getAsJsonObject().has("ClinVar_STATUS")
                            ? null
                            : jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_STATUS")))
                    .ClinVar_DBID(!wgsDataLocation.getAsJsonObject().has("ClinVar_DBID")
                            ? null
                            : jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("ClinVar_DBID")))
                    .ExAC_NFE(!wgsDataLocation.getAsJsonObject().has("ExAC_NFE")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_NFE").getAsString())
                    .ExAC_FIN(!wgsDataLocation.getAsJsonObject().has("ExAC_FIN")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_FIN").getAsString())
                    .G1000_ALL(!wgsDataLocation.getAsJsonObject().has("1000G_ALL")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_ALL").getAsString())
                    .G1000_SAS(!wgsDataLocation.getAsJsonObject().has("1000G_SAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_SAS").getAsString())
                    .G1000_EAS(!wgsDataLocation.getAsJsonObject().has("1000G_EAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_EAS").getAsString())
                    .G1000_AFR(!wgsDataLocation.getAsJsonObject().has("1000G_AFR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_AFR").getAsString())
                    .ExAC_SAS(!wgsDataLocation.getAsJsonObject().has("ExAC_SAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_SAS").getAsString())
                    .ExAC_EAS(!wgsDataLocation.getAsJsonObject().has("ExAC_EAS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_EAS").getAsString())
                    .ExAC_AMR(!wgsDataLocation.getAsJsonObject().has("ExAC_AMR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_AMR").getAsString())
                    .ExAC_AFR(!wgsDataLocation.getAsJsonObject().has("ExAC_AFR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_AFR").getAsString())
                    .ExAC_Freq(!wgsDataLocation.getAsJsonObject().has("ExAC_Freq")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ExAC_Freq").getAsString())
                    .End(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("End").getAsString())
                    .Start(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Start").getAsString())
                    .SiPhy_29way_logOdds(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("SiPhy_29way_logOdds").getAsString())
                    .FullAA(jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("FullAA")))
                    .Ref(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Ref").getAsString())
                    .GERP_RS(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GERP++_RS").getAsString())
                    .FATHMM(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("FATHMM").getAsString())
                    .NucleotideChange(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("NucleotideChange").getAsString())
                    .phyloP100way_vertebrate(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP100way_vertebrate").getAsString())
                    .Func(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP100way_vertebrate").getAsString())
                    .GWAS_PUBMED(!wgsDataLocation.getAsJsonObject().has("GWAS_PUBMED")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_PUBMED").getAsString())
                    .Transcript(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Transcript").getAsString())
                    .ESP6500si_AA(!wgsDataLocation.getAsJsonObject().has("ESP6500si_AA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ESP6500si_AA").getAsString())
                    .ESP6500si_EA(!wgsDataLocation.getAsJsonObject().has("ESP6500si_EA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("ESP6500si_EA").getAsString())
                    .G1000_EUR(!wgsDataLocation.getAsJsonObject().has("1000G_EUR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_EUR").getAsString())
                    .G1000_AMR(!wgsDataLocation.getAsJsonObject().has("1000G_AMR")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("1000G_AMR").getAsString())
                    .Chr_Start_Ref_Alt(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Chr_Start_Ref_Alt").getAsString())
                    .AA(!wgsDataLocation.getAsJsonObject().has("AA")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("AA").getAsString())
                    .PopFreqMax(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("PopFreqMax").getAsString())
                    .FATHMM_Pred(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("FATHMM_Pred").getAsString())
                    .wgRna(!wgsDataLocation.getAsJsonObject().has("wgRna")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("wgRna").getAsString())
                    .Gene(jsonArrayToStringList(wgsDataLocation.getAsJsonObject().getAsJsonArray("Gene")))
                    .phyloP46way_placental(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("phyloP46way_placental").getAsString())
                    .key(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("_key").getAsString())
                    .targetScanS(!wgsDataLocation.getAsJsonObject().has("targetScanS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("targetScanS").getAsString())
                    .Chr(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Chr").getAsString())
                    .COSMIC_ID(!wgsDataLocation.getAsJsonObject().has("COSMIC_ID")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .alt(wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("Alt").getAsString())
                    .GWAS_DIS(!wgsDataLocation.getAsJsonObject().has("GWAS_DIS")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_DIS").getAsString())
                    .GWAS_SNP(!wgsDataLocation.getAsJsonObject().has("GWAS_SNP")
                            ? null
                            : wgsDataLocation.getAsJsonObject().getAsJsonPrimitive("GWAS_SNP").getAsString())
                    .build());

        }
        return molecularMatchWGSadataLocationList;
    }

    @NotNull
    private static List<MolecularMatchGRch37Location> createGRCH37Location(@NotNull JsonArray arrayLocation) {
        List<MolecularMatchGRch37Location> gRch37LocationList = Lists.newArrayList();

        for (JsonElement location : arrayLocation) {
            Set<String> keysLocation = location.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES.contains(keysLocation.size())) {
                LOGGER.warn("Found " + keysLocation.size() + " in molecular match grch37 location rather than the expected "
                        + EXPECTED_MOLECULARMATCH_LOCATIONGRCH37_SIZES);
                LOGGER.warn(keysLocation);
            }

            gRch37LocationList.add(ImmutableMolecularMatchGRch37Location.builder()
                    .compositeKey(location.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .ref(location.getAsJsonObject().get("ref").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .stop(location.getAsJsonObject().get("stop").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(location.getAsJsonObject().get("start").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(location.getAsJsonObject().get("chr").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .alt(location.getAsJsonObject().get("alt").isJsonNull()
                            ? null
                            : location.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .validated(location.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcriptConsequences(createConsequencesGRCH37(location.getAsJsonObject().getAsJsonArray("transcript_consequences")))
                    .strand(location.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .build());
        }
        return gRch37LocationList;

    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequencesGRCH37> createConsequencesGRCH37(
            @NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequencesGRCH37> transcriptConsequencesGRCH37List = Lists.newArrayList();
        for (JsonElement transcriptConsequences : arrayTranscriptConsequence) {
            Set<String> keysTranscriptConsequences = transcriptConsequences.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES.contains(keysTranscriptConsequences.size())) {
                LOGGER.warn("Found " + keysTranscriptConsequences.size()
                        + " in molecular match transcript consequences grch37 rather than the expected "
                        + EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES__GRCH37_SIZES);
                LOGGER.warn(keysTranscriptConsequences);
            }

            transcriptConsequencesGRCH37List.add(ImmutableMolecularMatchTranscriptConsequencesGRCH37.builder()
                    .aminoAcidChange(transcriptConsequences.getAsJsonObject().get("amino_acid_change").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .txSites(jsonArrayToStringList(transcriptConsequences.getAsJsonObject().getAsJsonArray("txSites")))
                    .exonNumber(transcriptConsequences.getAsJsonObject().get("exonNumber").isJsonNull()
                            ? null
                            : createArrayExonNumber(transcriptConsequences))
                    .intronNumber(transcriptConsequences.getAsJsonObject().get("intronNumber").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .transcript(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(transcriptConsequences.getAsJsonObject().get("cdna").isJsonNull()
                            ? null
                            : transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return transcriptConsequencesGRCH37List;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequence> createTranscriptConsequence(@NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();

        for (JsonElement transcriptConsequence : arrayTranscriptConsequence) {
            Set<String> keysTranscriptConsequence = transcriptConsequence.getAsJsonObject().keySet();
            if (!EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES.contains(keysTranscriptConsequence.size())) {
                LOGGER.warn(
                        "Found " + keysTranscriptConsequence.size() + " in molecular match transcript consequence rather than the expected "
                                + EXPECTED_MOLECULARMATCH_TRANSCRIPT_CONSEQUENCES_SIZES);
                LOGGER.warn(keysTranscriptConsequence);
            }

            transcriptConsequenceList.add(ImmutableMolecularMatchTranscriptConsequence.builder()
                    .aminoAcidChange(
                            !transcriptConsequence.getAsJsonObject().has("amino_acid_change") || transcriptConsequence.getAsJsonObject()
                                    .get("amino_acid_change")
                                    .isJsonNull() ? null : transcriptConsequence.getAsJsonObject().get("amino_acid_change").getAsString())
                    .compositeKey(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .intronNumber(transcriptConsequence.getAsJsonObject().get("intronNumber").isJsonNull()
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(transcriptConsequence.getAsJsonObject().get("exonNumber").isJsonNull()
                            ? null
                            : createArrayExonNumber(transcriptConsequence))
                    .suppress(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .stop(!transcriptConsequence.getAsJsonObject().has("stop")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .custom(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .start(!transcriptConsequence.getAsJsonObject().has("start")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(!transcriptConsequence.getAsJsonObject().has("chr")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .validated(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcript(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(!transcriptConsequence.getAsJsonObject().has("cdna")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .referenceGenome(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(!transcriptConsequence.getAsJsonObject().has("ref")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .alt(!transcriptConsequence.getAsJsonObject().has("alt")
                            ? null
                            : transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .build());
        }
        return transcriptConsequenceList;
    }

    @NotNull
    private static JaxTrials createJaxTrials(@NotNull JsonObject objectJaxTrials) {
        return ImmutableJaxTrials.builder()
                .indications(createJaxTrialsIndications(objectJaxTrials.getAsJsonArray("indications")))
                .title(objectJaxTrials.getAsJsonPrimitive("title").getAsString())
                .gender(objectJaxTrials.get("gender").isJsonNull() ? null : objectJaxTrials.getAsJsonPrimitive("gender").getAsString())
                .nctId(objectJaxTrials.getAsJsonPrimitive("nctId").getAsString())
                .sponsors(objectJaxTrials.getAsJsonPrimitive("sponsors").getAsString())
                .recruitment(objectJaxTrials.getAsJsonPrimitive("recruitment").getAsString())
                .variantRequirements(objectJaxTrials.getAsJsonPrimitive("variantRequirements").getAsString())
                .updateDate(objectJaxTrials.getAsJsonPrimitive("updateDate").getAsString())
                .phase(objectJaxTrials.getAsJsonPrimitive("phase").getAsString())
                .variantRequirementDetails(createJaxTrialsVariantRequirementsDetails(objectJaxTrials.getAsJsonArray(
                        "variantRequirementDetails")))
                .therapies(createJaxTrialsTherapies(objectJaxTrials.getAsJsonArray("therapies")))
                .build();
    }

    @NotNull
    private static List<JaxTrialsIndications> createJaxTrialsIndications(@NotNull JsonArray arrayIndications) {
        List<JaxTrialsIndications> indicationsList = Lists.newArrayList();

        for (JsonElement indications : arrayIndications) {
            Set<String> keysIndications = indications.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES.contains(keysIndications.size())) {
                LOGGER.warn("Found " + keysIndications.size() + " in jax trials indications rather than the expected "
                        + EXPECTED_JAX_TRIALS_INDICATIONS_ELEMENT_SIZES);
                LOGGER.warn(keysIndications);
            }
            indicationsList.add(ImmutableJaxTrialsIndications.builder()
                    .source(indications.getAsJsonObject().getAsJsonPrimitive("source").getAsString())
                    .id(indications.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(indications.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return indicationsList;
    }

    @NotNull
    private static List<JaxTrialsVariantRequirementDetails> createJaxTrialsVariantRequirementsDetails(
            @NotNull JsonArray arrarVariantRequirementDetails) {
        List<JaxTrialsVariantRequirementDetails> variantRequirementDetailsList = Lists.newArrayList();
        for (JsonElement variantRequirementDetails : arrarVariantRequirementDetails) {
            Set<String> keysRequirementDetails = variantRequirementDetails.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES.contains(keysRequirementDetails.size())) {
                LOGGER.warn("Found " + keysRequirementDetails.size() + " in jax trials requirement details rather than the expected "
                        + EXPECTED_JAX_TRIALS_VARIANTREQUIREMENTDETAILS_ELEMENT_SIZES);
                LOGGER.warn(keysRequirementDetails);
            }
            variantRequirementDetailsList.add(ImmutableJaxTrialsVariantRequirementDetails.builder()
                    .molecularProfiles(createJaxTrialsMolecularProfile(variantRequirementDetails.getAsJsonObject()
                            .getAsJsonObject("molecularProfile")))
                    .requirementType(variantRequirementDetails.getAsJsonObject().getAsJsonPrimitive("requirementType").getAsString())
                    .build());
        }
        return variantRequirementDetailsList;
    }

    @NotNull
    private static List<JaxTrialsMolecularProfile> createJaxTrialsMolecularProfile(@NotNull JsonObject objectMolecularProfile) {
        List<JaxTrialsMolecularProfile> molecularProfileList = Lists.newArrayList();

        Set<String> keysMolecularProfile = objectMolecularProfile.keySet();
        if (!EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES.contains(keysMolecularProfile.size())) {
            LOGGER.warn("Found " + keysMolecularProfile.size() + " in jax trials molecular profile rather than the expected "
                    + EXPECTED_JAX_TRIALS_MOLECULAIRPROFILE_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularProfile);
        }

        molecularProfileList.add(ImmutableJaxTrialsMolecularProfile.builder()
                .profileName(objectMolecularProfile.getAsJsonPrimitive("profileName").getAsString())
                .id(objectMolecularProfile.getAsJsonPrimitive("id").getAsString())
                .build());
        return molecularProfileList;
    }

    @NotNull
    private static List<JaxTrialsTherapies> createJaxTrialsTherapies(@NotNull JsonArray arrayTherapies) {
        List<JaxTrialsTherapies> jaxTrialsTherapiesList = Lists.newArrayList();
        for (JsonElement therapies : arrayTherapies) {
            Set<String> keysTherapies = therapies.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES.contains(keysTherapies.size())) {
                LOGGER.warn("Found " + keysTherapies.size() + " in jax trials therapies rather than the expected "
                        + EXPECTED_JAX_TRIALS_THERAPIES_ELEMENT_SIZES);
                LOGGER.warn(keysTherapies);
            }

            jaxTrialsTherapiesList.add(ImmutableJaxTrialsTherapies.builder()
                    .id(therapies.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .therapyName(therapies.getAsJsonObject().getAsJsonPrimitive("therapyName").getAsString())
                    .build());
        }
        return jaxTrialsTherapiesList;
    }

    @NotNull
    private static Jax createJax(@NotNull JsonObject objectJax) {
        return ImmutableJax.builder()
                .responseType(objectJax.getAsJsonPrimitive("responseType").getAsString())
                .approvalStatus(objectJax.getAsJsonPrimitive("approvalStatus").getAsString())
                .molecularProfile(createMolecularProfile(objectJax.getAsJsonObject("molecularProfile")))
                .therapy(createJaxTherapy(objectJax.getAsJsonObject("therapy")))
                .evidenceType(objectJax.getAsJsonPrimitive("evidenceType").getAsString())
                .indications(createJaxIndications(objectJax.getAsJsonObject("indication")))
                .efficacyEvidence(objectJax.getAsJsonPrimitive("efficacyEvidence").getAsString())
                .references(objectJax.get("references").isJsonNull() ? null : createJaxReferences(objectJax.getAsJsonArray("references")))
                .id(objectJax.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static JaxMolecularProfile createMolecularProfile(@NotNull JsonObject objectMolecularProfile) {
        Set<String> keysMolecularProfile = objectMolecularProfile.keySet();
        if (!EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES.contains(keysMolecularProfile.size())) {
            LOGGER.warn("Found " + keysMolecularProfile.size() + " in jax molecular profile rather than the expected "
                    + EXPECTED_JAX_MOLECULAR_PROFILE_ELEMENT_SIZES);
            LOGGER.warn(keysMolecularProfile);
        }
        return ImmutableJaxMolecularProfile.builder()
                .profileName(objectMolecularProfile.getAsJsonPrimitive("profileName").getAsString())
                .id(objectMolecularProfile.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static JaxTherapy createJaxTherapy(@NotNull JsonObject objectTherapy) {
        Set<String> keysTherapy = objectTherapy.keySet();
        if (!EXPECTED_JAX_THERAPY_ELEMENT_SIZES.contains(keysTherapy.size())) {
            LOGGER.warn("Found " + keysTherapy.size() + " in jax therapy rather than the expected " + EXPECTED_JAX_THERAPY_ELEMENT_SIZES);
            LOGGER.warn(keysTherapy);
        }
        return ImmutableJaxTherapy.builder()
                .id(objectTherapy.getAsJsonPrimitive("id").getAsString())
                .therapyName(objectTherapy.getAsJsonPrimitive("therapyName").getAsString())
                .build();
    }

    @NotNull
    private static JaxIndications createJaxIndications(@NotNull JsonObject objectIndications) {
        Set<String> keysIndications = objectIndications.keySet();
        if (!EXPECTED_JAX_INDICATIONS_SIZES.contains(keysIndications.size())) {
            LOGGER.warn(
                    "Found " + keysIndications.size() + " in jax indications rather than the expected " + EXPECTED_JAX_INDICATIONS_SIZES);
            LOGGER.warn(keysIndications);
        }
        return ImmutableJaxIndications.builder()
                .source(objectIndications.getAsJsonPrimitive("source").getAsString())
                .id(objectIndications.getAsJsonPrimitive("id").getAsString())
                .name(objectIndications.getAsJsonPrimitive("name").getAsString())
                .build();
    }

    @NotNull
    private static List<JaxReferences> createJaxReferences(@NotNull JsonArray objectReferences) {
        List<JaxReferences> listReferences = Lists.newArrayList();
        for (JsonElement references : objectReferences) {
            Set<String> keysReferences = references.getAsJsonObject().keySet();
            if (!EXPECTED_JAX_REFERENCES_ELEMENT_SIZES.contains(keysReferences.size())) {
                LOGGER.warn("Found " + keysReferences.size() + " in jax references rather than the expected "
                        + EXPECTED_JAX_REFERENCES_ELEMENT_SIZES);
                LOGGER.warn(keysReferences);
            }
            listReferences.add(ImmutableJaxReferences.builder()
                    .url(references.getAsJsonObject().getAsJsonPrimitive("url").getAsString())
                    .id(references.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .pubMedId(references.getAsJsonObject().get("pubMedId").isJsonNull()
                            ? null
                            : references.getAsJsonObject().getAsJsonPrimitive("pubMedId").getAsString())
                    .title(references.getAsJsonObject().getAsJsonPrimitive("title").getAsString())
                    .build());
        }
        return listReferences;
    }

    @NotNull
    private static Oncokb createOncoKbBiological(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().oncoKbBiological(createBiologicalOncoKb(objectOncoKb.getAsJsonObject("biological"))).build();
    }

    @NotNull
    private static Oncokb createOncoKbClinical(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().oncoKbClinical(createClinicalOncoKb(objectOncoKb.getAsJsonObject("clinical"))).build();
    }

    @NotNull
    private static OncoKbClinical createClinicalOncoKb(@NotNull JsonObject objectClinical) {
        Set<String> keysClinical = objectClinical.keySet();
        if (!EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES.contains(keysClinical.size())) {
            LOGGER.warn("Found " + keysClinical.size() + " in oncokb clinical rather than the expected "
                    + EXPECTED_ONCOKB_CLINICAL_ELEMENT_SIZES);
            LOGGER.warn(keysClinical);
        }
        return ImmutableOncoKbClinical.builder()
                .RefSeq(objectClinical.getAsJsonPrimitive("RefSeq").getAsString())
                .level(objectClinical.getAsJsonPrimitive("level").getAsString())
                .Isoform(objectClinical.getAsJsonPrimitive("Isoform").getAsString())
                .oncokbVariant(createVariantOncoKb(objectClinical.getAsJsonObject("variant")))
                .entrezGeneID(objectClinical.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .drugPmids(objectClinical.getAsJsonPrimitive("drugPmids").getAsString())
                .cancerType(objectClinical.getAsJsonPrimitive("cancerType").getAsString())
                .drug(objectClinical.getAsJsonPrimitive("drug").getAsString())
                .gene(objectClinical.getAsJsonPrimitive("gene").getAsString())
                .levelLabel(objectClinical.getAsJsonPrimitive("level_label").getAsString())
                .oncoKbDrugAbstracts(createDrugsAbstracts(objectClinical.getAsJsonArray("drugAbstracts")))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstracts> createDrugsAbstracts(@NotNull JsonArray arrayDrugsAbstracts) {

        List<OncoKbDrugAbstracts> listDrugsabstracts = Lists.newArrayList();
        for (JsonElement drugAbstracts : arrayDrugsAbstracts) {
            Set<String> keysBiological = drugAbstracts.getAsJsonObject().keySet();

            if (!EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES.contains(keysBiological.size())) {
                LOGGER.warn("Found " + keysBiological.size() + " in oncokb drugs abstracts rather than the expected "
                        + EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES);
                LOGGER.warn(keysBiological);
            }
            listDrugsabstracts.add(ImmutableOncoKbDrugAbstracts.builder()
                    .text(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("text").getAsString())
                    .link(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .build());
        }
        return listDrugsabstracts;
    }

    @NotNull
    private static OncoKbBiological createBiologicalOncoKb(@NotNull JsonObject objectBiological) {
        Set<String> keysBiological = objectBiological.keySet();
        if (!EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES.contains(keysBiological.size())) {
            LOGGER.warn("Found " + keysBiological.size() + " in oncokb biological rather than the expected "
                    + EXPECTED_ONCOKB_BIOLOGICAL_ELEMENT_SIZES);
            LOGGER.warn(keysBiological);
        }

        return ImmutableOncoKbBiological.builder()
                .mutationEffectPmids(objectBiological.getAsJsonPrimitive("mutationEffectPmids").getAsString())
                .Isoform(objectBiological.getAsJsonPrimitive("Isoform").getAsString())
                .oncokbVariant(createVariantOncoKb(objectBiological.getAsJsonObject("variant")))
                .entrezGeneID(objectBiological.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .oncogenic(objectBiological.getAsJsonPrimitive("oncogenic").getAsString())
                .mutationEffect(objectBiological.getAsJsonPrimitive("mutationEffect").getAsString())
                .RefSeq(objectBiological.getAsJsonPrimitive("RefSeq").getAsString())
                .gene(objectBiological.getAsJsonPrimitive("gene").getAsString())
                .mutationEffectAbstracts(objectBiological.getAsJsonPrimitive("mutationEffectAbstracts").getAsString())
                .build();
    }

    @NotNull
    private static OncokbVariant createVariantOncoKb(@NotNull JsonObject objectVariant) {
        Set<String> keysVariant = objectVariant.keySet();

        if (!EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn(
                    "Found " + keysVariant.size() + " in oncokb variant rather than the expected" + EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }

        return ImmutableOncokbVariant.builder()
                .variantResidues(objectVariant.get("variantResidues").isJsonNull()
                        ? null
                        : objectVariant.getAsJsonPrimitive("variantResidues").getAsString())
                .proteinStart(objectVariant.getAsJsonPrimitive("proteinStart").getAsString())
                .name(objectVariant.getAsJsonPrimitive("name").getAsString())
                .proteinEnd(objectVariant.getAsJsonPrimitive("proteinEnd").getAsString())
                .refResidues(objectVariant.get("refResidues").isJsonNull()
                        ? null
                        : objectVariant.getAsJsonPrimitive("refResidues").getAsString())
                .alteration(objectVariant.getAsJsonPrimitive("alteration").getAsString())
                .oncoKbConsequence(createConsequenceOncokb(objectVariant.getAsJsonObject("consequence")))
                .oncokbGene(createGeneOncoKb(objectVariant.getAsJsonObject("gene")))
                .build();
    }

    @NotNull
    private static OncoKbConsequence createConsequenceOncokb(@NotNull JsonObject objectConsequence) {
        Set<String> keysConsequence = objectConsequence.keySet();

        if (!EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES.contains(keysConsequence.size())) {
            LOGGER.warn("Found " + keysConsequence.size() + " in oncokb consequence rather than the expected "
                    + EXPECTED_ONCOKB_CONSEQUENCE_ELEMENT_SIZES);
            LOGGER.warn(keysConsequence);
        }

        return ImmutableOncoKbConsequence.builder()
                .term(objectConsequence.getAsJsonPrimitive("term").getAsString())
                .description(objectConsequence.getAsJsonPrimitive("description").getAsString())
                .isGenerallyTruncating(objectConsequence.getAsJsonPrimitive("isGenerallyTruncating").getAsString())
                .build();
    }

    @NotNull
    private static OncokbGene createGeneOncoKb(@NotNull JsonObject objectGene) {
        Set<String> keysGene = objectGene.keySet();

        if (!EXPECTED_ONCOKB_GENE_ELEMENT_SIZES.contains(keysGene.size())) {
            LOGGER.warn("Found " + keysGene.size() + " in oncokb gene rather than the expected " + EXPECTED_ONCOKB_GENE_ELEMENT_SIZES);
            LOGGER.warn(keysGene);
        }

        return ImmutableOncokbGene.builder()
                .oncogene(objectGene.getAsJsonPrimitive("oncogene").getAsString())
                .name(objectGene.getAsJsonPrimitive("name").getAsString())
                .hugoSymbol(objectGene.getAsJsonPrimitive("hugoSymbol").getAsString())
                .curatedRefSeq(objectGene.get("curatedRefSeq").isJsonNull()
                        ? null
                        : objectGene.getAsJsonPrimitive("curatedRefSeq").getAsString())
                .entrezGeneId(objectGene.getAsJsonPrimitive("entrezGeneId").getAsString())
                .geneAliases(Lists.newArrayList(jsonArrayToStringList(objectGene.getAsJsonArray("geneAliases"))))
                .tsg(objectGene.getAsJsonPrimitive("tsg").getAsString())
                .curatedIsoform(objectGene.get("curatedIsoform").isJsonNull()
                        ? null
                        : objectGene.getAsJsonPrimitive("curatedIsoform").getAsString())
                .build();
    }

    @NotNull
    private static Pmkb createPmkb(@NotNull JsonObject objectPmkb) {
        JsonObject tumor = objectPmkb.getAsJsonObject("tumor");
        JsonArray tissue = objectPmkb.getAsJsonArray("tissues");
        JsonObject variant = objectPmkb.getAsJsonObject("variant");

        return ImmutablePmkb.builder().tumor(createTumor(tumor)).tissue(createTissue(tissue)).variant(createVariantPmkb(variant)).build();
    }

    @NotNull
    private static List<PmkbTumor> createTumor(@NotNull JsonObject tumor) {
        Set<String> keysTumor = tumor.keySet();
        if (!EXPECTED_PMKB_TUMOR_ELEMENT_SIZES.contains(keysTumor.size())) {
            LOGGER.warn("Found " + keysTumor.size() + " in pmkb tumor rather than the expected " + EXPECTED_PMKB_TUMOR_ELEMENT_SIZES);
            LOGGER.warn(keysTumor);
        }

        List<PmkbTumor> listTumor = Lists.newArrayList();
        listTumor.add(ImmutablePmkbTumor.builder()
                .id(tumor.getAsJsonPrimitive("id").getAsString())
                .name(tumor.getAsJsonPrimitive("name").getAsString())
                .build());

        return listTumor;
    }

    @NotNull
    private static List<PmkbTissue> createTissue(@NotNull JsonArray tissues) {

        List<PmkbTissue> listTissue = Lists.newArrayList();
        for (JsonElement tissue : tissues) {
            Set<String> keysTissue = tissue.getAsJsonObject().keySet();
            if (!EXPECTED_PMKB_TISSUE_ELEMENT_SIZES.contains(keysTissue.size())) {
                LOGGER.warn(
                        "Found " + keysTissue.size() + " in pmkb tissue rather than the expected " + EXPECTED_PMKB_TISSUE_ELEMENT_SIZES);
                LOGGER.warn(keysTissue);
            }
            listTissue.add(ImmutablePmkbTissue.builder()
                    .id(tissue.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(tissue.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return listTissue;
    }

    @NotNull
    private static List<PmkbVariant> createVariantPmkb(@NotNull JsonObject variant) {
        Set<String> keysVariant = variant.keySet();
        if (!EXPECTED_PMKB_VARIANT_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn("Found " + keysVariant.size() + " in pmkb variant rather than the expected " + EXPECTED_PMKB_VARIANT_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }
        return Lists.newArrayList(ImmutablePmkbVariant.builder()
                .aminoAcidChange(variant.get("amino_acid_change").isJsonNull() ? null : variant.get("amino_acid_change").getAsString())
                .germline(variant.get("germline").isJsonNull() ? null : variant.get("germline").getAsString())
                .partnerGene(variant.get("partner_gene").isJsonNull() ? null : variant.get("partner_gene").getAsString())
                .codons(variant.get("codons").isJsonNull() ? null : variant.get("codons").getAsString())
                .description(variant.get("description").isJsonNull() ? null : variant.get("description").getAsString())
                .exons(variant.get("exons").isJsonNull() ? null : variant.get("exons").getAsString())
                .notes(variant.get("notes").isJsonNull() ? null : variant.get("notes").getAsString())
                .cosmic(variant.get("cosmic").isJsonNull() ? null : variant.get("cosmic").getAsString())
                .effect(variant.get("effect").isJsonNull() ? null : variant.get("effect").getAsString())
                .cnvType(variant.get("cnv_type").isJsonNull() ? null : variant.get("cnv_type").getAsString())
                .id(variant.get("id").isJsonNull() ? null : variant.get("id").getAsString())
                .cytoband(variant.get("cytoband").isJsonNull() ? null : variant.get("cytoband").getAsString())
                .variantType(variant.get("variant_type").isJsonNull() ? null : variant.get("variant_type").getAsString())
                .dnaChange(variant.get("dna_change").isJsonNull() ? null : variant.get("dna_change").getAsString())
                .coordinates(variant.get("coordinates").isJsonNull() ? null : variant.get("coordinates").getAsString())
                .chromosomeBasedCnv(variant.get("chromosome_based_cnv").isJsonNull()
                        ? null
                        : variant.get("chromosome_based_cnv").getAsString())
                .gene(createGene(variant))
                .transcript(variant.getAsJsonPrimitive("transcript").getAsString())
                .descriptionType(variant.get("description_type").isJsonNull() ? null : variant.get("description_type").getAsString())
                .chromosome(variant.get("chromosome").isJsonNull() ? null : variant.get("chromosome").getAsString())
                .name(variant.get("name").isJsonNull() ? null : variant.get("name").getAsString())
                .build());
    }

    @NotNull
    private static List<PmkbGene> createGene(@NotNull JsonObject variant) {
        JsonObject gene = variant.getAsJsonObject("gene");

        Set<String> keysgene = gene.keySet();
        if (!EXPECTED_PMKB_GENE_ELEMENT_SIZES.contains(keysgene.size())) {
            LOGGER.warn("Found " + keysgene.size() + " in pmkb gene rather than the expected " + EXPECTED_PMKB_GENE_ELEMENT_SIZES);
            LOGGER.warn(keysgene);
        }

        List<PmkbGene> listGene = Lists.newArrayList();
        listGene.add(ImmutablePmkbGene.builder()
                .description(gene.get("description").isJsonNull() ? null : gene.getAsJsonPrimitive("description").getAsString())
                .createdAt(gene.getAsJsonPrimitive("created_at").getAsString())
                .updatedAt(gene.getAsJsonPrimitive("updated_at").getAsString())
                .activeInd(gene.getAsJsonPrimitive("active_ind").getAsString())
                .externalId(gene.getAsJsonPrimitive("external_id").getAsString())
                .id(gene.getAsJsonPrimitive("id").getAsString())
                .name(gene.getAsJsonPrimitive("name").getAsString())
                .build());
        return listGene;
    }

    @NotNull
    private static Sage createSage(@NotNull JsonObject objectSage) {
        return ImmutableSage.builder()
                .entrezId(objectSage.getAsJsonPrimitive("entrez_id").getAsString())
                .clinicalManifestation(objectSage.getAsJsonPrimitive("clinical_manifestation").getAsString())
                .publicationUrl(objectSage.getAsJsonPrimitive("publication_url").getAsString())
                .germlineOrSomatic(objectSage.getAsJsonPrimitive("germline_or_somatic").getAsString())
                .evidenceLabel(objectSage.getAsJsonPrimitive("evidence_label").getAsString())
                .drugLabels(objectSage.getAsJsonPrimitive("drug_labels").getAsString())
                .responseType(objectSage.getAsJsonPrimitive("response_type").getAsString())
                .gene(objectSage.getAsJsonPrimitive("gene").getAsString())
                .build();
    }

    @NotNull
    private static BRCA createBRCA(@NotNull JsonObject objectBrca) {
        return ImmutableBRCA.builder().brcApart1(createBRCAPart1(objectBrca)).brcApart2(createBRCAPart2(objectBrca)).build();
    }

    @NotNull
    private static BRCApart1 createBRCAPart1(@NotNull JsonObject objectBrca) {
        return ImmutableBRCApart1.builder()
                .Variant_frequency_LOVD(objectBrca.getAsJsonPrimitive("Variant_frequency_LOVD").getAsString())
                .Allele_frequency_FIN_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_FIN_ExAC").getAsString())
                .ClinVarAccession_ENIGMA(objectBrca.getAsJsonPrimitive("ClinVarAccession_ENIGMA").getAsString())
                .Homozygous_count_AFR_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_AFR_ExAC").getAsString())
                .BX_ID_ExAC(objectBrca.getAsJsonPrimitive("BX_ID_ExAC").getAsString())
                .Variant_in_LOVD(objectBrca.getAsJsonPrimitive("Variant_in_LOVD").getAsString())
                .Allele_frequency_AFR_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_AFR_ExAC").getAsString())
                .Chr(objectBrca.getAsJsonPrimitive("Chr").getAsString())
                .BX_ID_ENIGMA(objectBrca.getAsJsonPrimitive("BX_ID_ENIGMA").getAsString())
                .Co_occurrence_LR_exLOVD(objectBrca.getAsJsonPrimitive("Co_occurrence_LR_exLOVD").getAsString())
                .Homozygous_count_EAS_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_EAS_ExAC").getAsString())
                .Submitter_ClinVar(objectBrca.getAsJsonPrimitive("Submitter_ClinVar").getAsString())
                .Allele_frequency_EAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_EAS_ExAC").getAsString())
                .Hg37_End(objectBrca.getAsJsonPrimitive("Hg37_End").getAsString())
                .Submitters_LOVD(objectBrca.getAsJsonPrimitive("Submitters_LOVD").getAsString())
                .Clinical_classification_BIC(objectBrca.getAsJsonPrimitive("Clinical_classification_BIC").getAsString())
                .Homozygous_count_NFE_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_NFE_ExAC").getAsString())
                .Allele_count_SAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_SAS_ExAC").getAsString())
                .Method_ClinVar(objectBrca.getAsJsonPrimitive("Method_ClinVar").getAsString())
                .Allele_count_NFE_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_NFE_ExAC").getAsString())
                .Pathogenicity_all(objectBrca.getAsJsonPrimitive("Pathogenicity_all").getAsString())
                .Germline_or_Somatic_BIC(objectBrca.getAsJsonPrimitive("Germline_or_Somatic_BIC").getAsString())
                .Homozygous_count_SAS_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_SAS_ExAC").getAsString())
                .BIC_Nomenclature(objectBrca.getAsJsonPrimitive("BIC_Nomenclature").getAsString())
                .Assertion_method_ENIGMA(objectBrca.getAsJsonPrimitive("Assertion_method_ENIGMA").getAsString())
                .Literature_source_exLOVD(objectBrca.getAsJsonPrimitive("Literature_source_exLOVD").getAsString())
                .Change_Type_id(objectBrca.getAsJsonPrimitive("Change_Type_id").getAsString())
                .Collection_method_ENIGMA(objectBrca.getAsJsonPrimitive("Collection_method_ENIGMA").getAsString())
                .Sum_family_LR_exLOVD(objectBrca.getAsJsonPrimitive("Sum_family_LR_exLOVD").getAsString())
                .HGVS_cDNA_LOVD(objectBrca.getAsJsonPrimitive("HGVS_cDNA_LOVD").getAsString())
                .Homozygous_count_FIN_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_FIN_ExAC").getAsString())
                .EAS_Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("EAS_Allele_frequency_1000_Genomes").getAsString())
                .Ethnicity_BIC(objectBrca.getAsJsonPrimitive("Ethnicity_BIC").getAsString())
                .Individuals_LOVD(objectBrca.getAsJsonPrimitive("Individuals_LOVD").getAsString())
                .Variant_in_ExAC(objectBrca.getAsJsonPrimitive("Variant_in_ExAC").getAsString())
                .URL_ENIGMA(objectBrca.getAsJsonPrimitive("URL_ENIGMA").getAsString())
                .Allele_Origin_ClinVar(objectBrca.getAsJsonPrimitive("Allele_Origin_ClinVar").getAsString())
                .Allele_frequency_AMR_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_AMR_ExAC").getAsString())
                .Variant_in_1000_Genomes(objectBrca.getAsJsonPrimitive("Variant_in_1000_Genomes").getAsString())
                .AFR_Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("AFR_Allele_frequency_1000_Genomes").getAsString())
                .BX_ID_exLOVD(objectBrca.getAsJsonPrimitive("BX_ID_exLOVD").getAsString())
                .Source(objectBrca.getAsJsonPrimitive("Source").getAsString())
                .Condition_ID_value_ENIGMA(objectBrca.getAsJsonPrimitive("Condition_ID_value_ENIGMA").getAsString())
                .HGVS_Protein(objectBrca.getAsJsonPrimitive("HGVS_Protein").getAsString())
                .Ref(objectBrca.getAsJsonPrimitive("Ref").getAsString())
                .Allele_number_AFR_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_AFR_ExAC").getAsString())
                .Allele_count_AFR_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_AFR_ExAC").getAsString())
                .BX_ID_LOVD(objectBrca.getAsJsonPrimitive("BX_ID_LOVD").getAsString())
                .Synonyms(objectBrca.getAsJsonPrimitive("Synonyms").getAsString())
                .Gene_Symbol(objectBrca.getAsJsonPrimitive("Gene_Symbol").getAsString())
                .Comment_on_clinical_significance_ENIGMA(objectBrca.getAsJsonPrimitive("Comment_on_clinical_significance_ENIGMA")
                        .getAsString())
                .Missense_analysis_prior_probability_exLOVD(objectBrca.getAsJsonPrimitive("Missense_analysis_prior_probability_exLOVD")
                        .getAsString())
                .Allele_number_FIN_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_FIN_ExAC").getAsString())
                .Posterior_probability_exLOVD(objectBrca.getAsJsonPrimitive("Posterior_probability_exLOVD").getAsString())
                .Polyphen_Score(objectBrca.getAsJsonPrimitive("Polyphen_Score").getAsString())
                .Reference_Sequence(objectBrca.getAsJsonPrimitive("Reference_Sequence").getAsString())
                .Allele_count_EAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_EAS_ExAC").getAsString())
                .Hg38_End(objectBrca.getAsJsonPrimitive("Hg38_End").getAsString())
                .HGVS_cDNA(objectBrca.getAsJsonPrimitive("HGVS_cDNA").getAsString())
                .Functional_analysis_technique_LOVD(objectBrca.getAsJsonPrimitive("Functional_analysis_technique_LOVD").getAsString())
                .SAS_Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("SAS_Allele_frequency_1000_Genomes").getAsString())
                .RNA_LOVD(objectBrca.getAsJsonPrimitive("RNA_LOVD").getAsString())
                .Combined_prior_probablility_exLOVD(objectBrca.getAsJsonPrimitive("Combined_prior_probablility_exLOVD").getAsString())
                .BX_ID_ClinVar(objectBrca.getAsJsonPrimitive("BX_ID_ClinVar").getAsString())
                .IARC_class_exLOVD(objectBrca.getAsJsonPrimitive("IARC_class_exLOVD").getAsString())
                .BX_ID_BIC(objectBrca.getAsJsonPrimitive("BX_ID_BIC").getAsString())
                .Sift_Prediction(objectBrca.getAsJsonPrimitive("Sift_Prediction").getAsString())
                .Allele_number_NFE_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_NFE_ExAC").getAsString())
                .Allele_origin_ENIGMA(objectBrca.getAsJsonPrimitive("Allele_origin_ENIGMA").getAsString())
                .Allele_number_OTH_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_OTH_ExAC").getAsString())
                .Hg36_End(objectBrca.getAsJsonPrimitive("Hg36_End").getAsString())
                .Allele_frequency_SAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_SAS_ExAC").getAsString())
                .Date_Last_Updated_ClinVar(objectBrca.getAsJsonPrimitive("Date_Last_Updated_ClinVar").getAsString())
                .Allele_number_EAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_EAS_ExAC").getAsString())
                .Allele_frequency_OTH_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_OTH_ExAC").getAsString())
                .Source_URL(objectBrca.getAsJsonPrimitive("Source_URL").getAsString())
                .SCV_ClinVar(objectBrca.getAsJsonPrimitive("SCV_ClinVar").getAsString())
                .Pathogenicity_expert(objectBrca.getAsJsonPrimitive("Pathogenicity_expert").getAsString())
                .Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("Allele_frequency_1000_Genomes").getAsString())
                .Functional_analysis_result_LOVD(objectBrca.getAsJsonPrimitive("Functional_analysis_result_LOVD").getAsString())
                .AMR_Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("AMR_Allele_frequency_1000_Genomes").getAsString())
                .Variant_in_ESP(objectBrca.getAsJsonPrimitive("Variant_in_ESP").getAsString())
                .Variant_in_BIC(objectBrca.getAsJsonPrimitive("Variant_in_BIC").getAsString())
                .Clinical_significance_ENIGMA(objectBrca.getAsJsonPrimitive("Clinical_significance_ENIGMA").getAsString())
                .Max_Allele_Frequency(objectBrca.getAsJsonPrimitive("Max_Allele_Frequency").getAsString())
                .Allele_count_AMR_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_AMR_ExAC").getAsString())
                .Variant_in_ENIGMA(objectBrca.getAsJsonPrimitive("Variant_in_ENIGMA").getAsString())
                .BX_ID_ESP(objectBrca.getAsJsonPrimitive("BX_ID_ESP").getAsString())
                .Patient_nationality_BIC(objectBrca.getAsJsonPrimitive("Patient_nationality_BIC").getAsString())
                .BX_ID_1000_Genomes(objectBrca.getAsJsonPrimitive("BX_ID_1000_Genomes").getAsString())
                .Genomic_Coordinate_hg37(objectBrca.getAsJsonPrimitive("Genomic_Coordinate_hg37").getAsString())
                .Genomic_Coordinate_hg36(objectBrca.getAsJsonPrimitive("Genomic_Coordinate_hg36").getAsString())
                .EUR_Allele_frequency_1000_Genomes(objectBrca.getAsJsonPrimitive("EUR_Allele_frequency_1000_Genomes").getAsString())
                .Number_of_family_member_carrying_mutation_BIC(objectBrca.getAsJsonPrimitive("Number_of_family_member_carrying_mutation_BIC")
                        .getAsString())
                .Segregation_LR_exLOVD(objectBrca.getAsJsonPrimitive("Segregation_LR_exLOVD").getAsString())
                .Allele_Frequency(objectBrca.getAsJsonPrimitive("Allele_Frequency").getAsString())
                .Minor_allele_frequency_percent_ESP(objectBrca.getAsJsonPrimitive("Minor_allele_frequency_percent_ESP").getAsString())
                .Allele_frequency_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_ExAC").getAsString())
                .Mutation_type_BIC(objectBrca.getAsJsonPrimitive("Mutation_type_BIC").getAsString())
                .Assertion_method_citation_ENIGMA(objectBrca.getAsJsonPrimitive("Assertion_method_citation_ENIGMA").getAsString())
                .Condition_ID_type_ENIGMA(objectBrca.getAsJsonPrimitive("Condition_ID_type_ENIGMA").getAsString())
                .Allele_count_OTH_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_OTH_ExAC").getAsString())
                .HGVS_protein_LOVD(objectBrca.getAsJsonPrimitive("HGVS_protein_LOVD").getAsString())
                .Variant_in_ClinVar(objectBrca.getAsJsonPrimitive("Variant_in_ClinVar").getAsString())
                .Clinical_importance_BIC(objectBrca.getAsJsonPrimitive("Clinical_importance_BIC").getAsString())
                .Discordant(objectBrca.getAsJsonPrimitive("Discordant").getAsString())
                .build();
    }

    @NotNull
    private static BRCApart2 createBRCAPart2(@NotNull JsonObject objectBrca) {
        return ImmutableBRCApart2.builder()
                .Allele_count_FIN_ExAC(objectBrca.getAsJsonPrimitive("Allele_count_FIN_ExAC").getAsString())
                .Condition_category_ENIGMA(objectBrca.getAsJsonPrimitive("Condition_category_ENIGMA").getAsString())
                .Allele_Frequency_ESP(objectBrca.getAsJsonPrimitive("Allele_Frequency_ESP").getAsString())
                .Homozygous_count_OTH_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_OTH_ExAC").getAsString())
                .Genetic_origin_LOVD(objectBrca.getAsJsonPrimitive("Genetic_origin_LOVD").getAsString())
                .id(objectBrca.getAsJsonPrimitive("id").getAsString())
                .Homozygous_count_AMR_ExAC(objectBrca.getAsJsonPrimitive("Homozygous_count_AMR_ExAC").getAsString())
                .Clinical_Significance_ClinVar(objectBrca.getAsJsonPrimitive("Clinical_Significance_ClinVar").getAsString())
                .AA_Allele_Frequency_ESP(objectBrca.getAsJsonPrimitive("AA_Allele_Frequency_ESP").getAsString())
                .Protein_Change(objectBrca.getAsJsonPrimitive("Protein_Change").getAsString())
                .Variant_in_exLOVD(objectBrca.getAsJsonPrimitive("Variant_in_exLOVD").getAsString())
                .EA_Allele_Frequency_ESP(objectBrca.getAsJsonPrimitive("EA_Allele_Frequency_ESP").getAsString())
                .HGVS_RNA(objectBrca.getAsJsonPrimitive("HGVS_RNA").getAsString())
                .Clinical_significance_citations_ENIGMA(objectBrca.getAsJsonPrimitive("Clinical_significance_citations_ENIGMA")
                        .getAsString())
                .Variant_effect_LOVD(objectBrca.getAsJsonPrimitive("Variant_effect_LOVD").getAsString())
                .Polyphen_Prediction(objectBrca.getAsJsonPrimitive("Polyphen_Prediction").getAsString())
                .Data_Release_id(objectBrca.getAsJsonPrimitive("Data_Release_id").getAsString())
                .Hg37_Start(objectBrca.getAsJsonPrimitive("Hg37_Start").getAsString())
                .Hg36_Start(objectBrca.getAsJsonPrimitive("Hg36_Start").getAsString())
                .Sift_Score(objectBrca.getAsJsonPrimitive("Sift_Score").getAsString())
                .Genomic_Coordinate_hg38(objectBrca.getAsJsonPrimitive("Genomic_Coordinate_hg38").getAsString())
                .Alt(objectBrca.getAsJsonPrimitive("Alt").getAsString())
                .Literature_citation_BIC(objectBrca.getAsJsonPrimitive("Literature_citation_BIC").getAsString())
                .Variant_haplotype_LOVD(objectBrca.getAsJsonPrimitive("Variant_haplotype_LOVD").getAsString())
                .Allele_frequency_NFE_ExAC(objectBrca.getAsJsonPrimitive("Allele_frequency_NFE_ExAC").getAsString())
                .Hg38_Start(objectBrca.getAsJsonPrimitive("Hg38_Start").getAsString())
                .Pos(objectBrca.getAsJsonPrimitive("Pos").getAsString())
                .Date_last_evaluated_ENIGMA(objectBrca.getAsJsonPrimitive("Date_last_evaluated_ENIGMA").getAsString())
                .Allele_number_SAS_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_SAS_ExAC").getAsString())
                .Allele_number_AMR_ExAC(objectBrca.getAsJsonPrimitive("Allele_number_AMR_ExAC").getAsString())
                .DBID_LOVD(objectBrca.getAsJsonPrimitive("DBID_LOVD").getAsString())
                .build();
    }

    @NotNull
    private static Cgi createCgi(@NotNull JsonObject objectCgi) {
        return ImmutableCgi.builder()
                .targeting(objectCgi.getAsJsonPrimitive("Targeting").getAsString())
                .source(objectCgi.getAsJsonPrimitive("Source").getAsString())
                .cDNA(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("cDNA"))))
                .primary_tumor_type(objectCgi.getAsJsonPrimitive("Primary Tumor type").getAsString())
                .individual_mutation(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("individual_mutation"))))
                .drugsFullName(objectCgi.getAsJsonPrimitive("Drug full name").getAsString())
                .curator(objectCgi.getAsJsonPrimitive("Curator").getAsString())
                .drug_family(objectCgi.getAsJsonPrimitive("Drug family").getAsString())
                .alteration(objectCgi.getAsJsonPrimitive("Alteration").getAsString())
                .drug(objectCgi.getAsJsonPrimitive("Drug").getAsString())
                .biomarker(objectCgi.getAsJsonPrimitive("Biomarker").getAsString())
                .gDNA(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("gDNA"))))
                .drug_status(objectCgi.getAsJsonPrimitive("Drug status").getAsString())
                .gene(objectCgi.getAsJsonPrimitive("Gene").getAsString())
                .transcript(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("transcript"))))
                .strand(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("strand"))))
                .info(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("info"))))
                .assay_type(objectCgi.getAsJsonPrimitive("Assay type").getAsString())
                .alteration_type(objectCgi.getAsJsonPrimitive("Alteration type").getAsString())
                .region(Lists.newArrayList(jsonArrayToStringList(objectCgi.getAsJsonArray("region"))))
                .evidence_level(objectCgi.getAsJsonPrimitive("Evidence level").getAsString())
                .association(objectCgi.getAsJsonPrimitive("Association").getAsString())
                .metastatic_Tumor_Type(objectCgi.getAsJsonPrimitive("Metastatic Tumor Type").getAsString())
                .build();
    }

    @NotNull
    private static List<Feature> createFeatures(@NotNull JsonObject viccEntryObject) {

        JsonArray arrayFeatures = viccEntryObject.getAsJsonArray("features");
        List<Feature> featureList = Lists.newArrayList();

        for (JsonElement elementFeature : arrayFeatures) {
            JsonObject objectFeatures = elementFeature.getAsJsonObject();
            Set<String> keysFeatures = objectFeatures.keySet();
            if (!EXPECTED_FEATURES_ELEMENT_SIZES.contains(keysFeatures.size())) {
                LOGGER.warn("Found " + keysFeatures.size() + " in feature rather than the expected " + EXPECTED_FEATURES_ELEMENT_SIZES);
                LOGGER.warn(keysFeatures);
            }
            featureList.add(ImmutableFeature.builder()
                    .name(objectFeatures.has("name") ? objectFeatures.getAsJsonPrimitive("name").getAsString() : null)
                    .biomarkerType(objectFeatures.has("biomarker_type")
                            ? objectFeatures.getAsJsonPrimitive("biomarker_type").getAsString()
                            : null)
                    .referenceName(objectFeatures.has("referenceName")
                            ? objectFeatures.getAsJsonPrimitive("referenceName").getAsString()
                            : null)
                    .chromosome(objectFeatures.has("chromosome") ? objectFeatures.getAsJsonPrimitive("chromosome").getAsString() : null)
                    .start(objectFeatures.has("start") && !objectFeatures.get("start").isJsonNull() ? objectFeatures.getAsJsonPrimitive(
                            "start").getAsString() : null)
                    .end(objectFeatures.has("end") && !objectFeatures.get("end").isJsonNull() ? objectFeatures.getAsJsonPrimitive("end")
                            .getAsString() : null)
                    .ref(objectFeatures.has("ref") && !objectFeatures.get("ref").isJsonNull() ? objectFeatures.getAsJsonPrimitive("ref")
                            .getAsString() : null)
                    .alt(objectFeatures.has("alt") && !objectFeatures.get("alt").isJsonNull() ? objectFeatures.getAsJsonPrimitive("alt")
                            .getAsString() : null)
                    .provenance(Lists.newArrayList())
                    .provenanceRule(objectFeatures.has("provenance_rule") ? objectFeatures.getAsJsonPrimitive("provenance_rule")
                            .getAsString() : null)
                    .geneSymbol(objectFeatures.has("geneSymbol") && !objectFeatures.get("geneSymbol").isJsonNull()
                            ? objectFeatures.getAsJsonPrimitive("geneSymbol").getAsString()
                            : null)
                    .synonyms(objectFeatures.has("synonyms") ? jsonArrayToStringList(objectFeatures.getAsJsonArray("synonyms")) : null)
                    .entrezId(objectFeatures.has("entrez_id") ? objectFeatures.getAsJsonPrimitive("entrez_id").getAsString() : null)
                    .sequenceOntology(objectFeatures.has("sequence_ontology") ? createSequenceOntology(objectFeatures.getAsJsonObject(
                            "sequence_ontology")) : null)
                    .links(objectFeatures.has("links") ? jsonArrayToStringList(objectFeatures.getAsJsonArray("links")) : null)
                    .description(objectFeatures.has("description") ? objectFeatures.getAsJsonPrimitive("description").getAsString() : null)
                    .build());
        }

        return featureList;
    }

    @NotNull
    private static SequenceOntology createSequenceOntology(@NotNull JsonObject objectSequenceOntology) {
        Set<String> keysSequenceOntology = objectSequenceOntology.keySet();
        if (!EXPECTED_SEQUENCE_ONTOLOGY_ELEMENT_SIZES.contains(keysSequenceOntology.size())) {
            LOGGER.warn("Found " + keysSequenceOntology.size() + " in sequence ontology rather than the expected "
                    + EXPECTED_SEQUENCE_ONTOLOGY_ELEMENT_SIZES);
            LOGGER.warn(keysSequenceOntology);
        }

        return ImmutableSequenceOntology.builder()
                .hierarchy(objectSequenceOntology.has("hierarchy")
                        ? jsonArrayToStringList(objectSequenceOntology.getAsJsonArray("hierarchy"))
                        : null)
                .soid(objectSequenceOntology.getAsJsonPrimitive("soid").getAsString())
                .parentSoid(objectSequenceOntology.getAsJsonPrimitive("parent_soid").getAsString())
                .name(objectSequenceOntology.getAsJsonPrimitive("name").getAsString())
                .parentName(objectSequenceOntology.getAsJsonPrimitive("parent_name").getAsString())
                .build();
    }

    @NotNull
    private static List<GeneIdentifier> createGeneIdentifiers(@NotNull JsonObject viccEntryObject) {
        JsonArray geneIdentifiers = viccEntryObject.getAsJsonArray("gene_identifiers");
        List<GeneIdentifier> listGeneIdentifiers = Lists.newArrayList();

        for (JsonElement elementGeneIdentifier : geneIdentifiers) {
            Set<String> keysGeneIdentifier = elementGeneIdentifier.getAsJsonObject().keySet();
            if (!EXPECTED_GENE_IDENTIFIERS_ELEMENT_SIZES.contains(keysGeneIdentifier.size())) {
                LOGGER.warn("Found " + keysGeneIdentifier.size() + " in gene identifier rather than the expected "
                        + EXPECTED_GENE_IDENTIFIERS_ELEMENT_SIZES);
                LOGGER.warn(keysGeneIdentifier);
            }
            listGeneIdentifiers.add(toGeneIdentifier(elementGeneIdentifier.getAsJsonObject()));
        }
        return listGeneIdentifiers;
    }

    @NotNull
    private static GeneIdentifier toGeneIdentifier(@NotNull JsonObject geneIdentifierObject) {
        return ImmutableGeneIdentifier.builder()
                .symbol(geneIdentifierObject.getAsJsonPrimitive("symbol").getAsString())
                .entrezId(geneIdentifierObject.getAsJsonPrimitive("entrez_id").getAsString())
                .ensemblGeneId(!geneIdentifierObject.get("ensembl_gene_id").isJsonNull() ? geneIdentifierObject.getAsJsonPrimitive(
                        "ensembl_gene_id").getAsString() : null)
                .build();
    }

    @NotNull
    private static Association createAssociation(@NotNull JsonObject associationObject) {
        return ImmutableAssociation.builder()
                .variantName(associationObject.has("variant_name") && associationObject.get("variant_name").isJsonArray()
                        ? Strings.EMPTY
                        : associationObject.has("variant_name") && associationObject.get("variant_name").isJsonPrimitive()
                                ? associationObject.getAsJsonPrimitive("variant_name").getAsString()
                                : null)
                .evidence(createEvidence(associationObject.getAsJsonArray("evidence")))
                .evidenceLevel(associationObject.has("evidence_level") ? associationObject.getAsJsonPrimitive("evidence_level")
                        .getAsString() : null)
                .evidenceLabel(
                        associationObject.has("evidence_label") && !associationObject.get("evidence_label").isJsonNull() ? associationObject
                                .getAsJsonPrimitive("evidence_label")
                                .getAsString() : null)
                .responseType(associationObject.has("response_type") && !associationObject.get("response_type").isJsonNull()
                        ? associationObject.getAsJsonPrimitive("response_type").getAsString()
                        : null)
                .drugLabels(associationObject.has("drug_labels") ? associationObject.getAsJsonPrimitive("drug_labels").getAsString() : null)
                .sourceLink(associationObject.has("source_link") ? associationObject.getAsJsonPrimitive("source_link").getAsString() : null)
                .publicationUrls(associationObject.has("publication_url") && associationObject.get("publication_url").isJsonPrimitive()
                        ? Lists.newArrayList(associationObject.getAsJsonPrimitive("publication_url").getAsString())
                        : associationObject.has("publication_url") && associationObject.get("publication_url").isJsonArray()
                                ? Lists.newArrayList(associationObject.getAsJsonArray("publication_url").getAsString())
                                : null)
                .phenotype(associationObject.has("phenotype") ? createPhenotype(associationObject.getAsJsonObject("phenotype")) : null)
                .description(associationObject.getAsJsonPrimitive("description").getAsString())
                .environmentalContexts(associationObject.get("environmentalContexts") != null
                        ? createEnvironmentalContexts(associationObject.getAsJsonArray("environmentalContexts"))
                        : null)
                .oncogenic(associationObject.has("oncogenic") ? associationObject.getAsJsonPrimitive("oncogenic").getAsString() : null)
                .build();
    }

    @NotNull
    private static List<EnvironmentalContext> createEnvironmentalContexts(@NotNull JsonArray arrayEnvironmentalContexts) {
        List<EnvironmentalContext> environmentalContexts = Lists.newArrayList();

        for (JsonElement elementEnvironmentContext : arrayEnvironmentalContexts) {
            JsonObject environmentContextObject = elementEnvironmentContext.getAsJsonObject();

            List<String> approvedCountries = Lists.newArrayList();
            if (environmentContextObject.has("approved_countries")) {
                for (JsonElement approvedCountriesElement : environmentContextObject.getAsJsonArray("approved_countries")) {
                    approvedCountries.add(approvedCountriesElement.getAsString());
                }
            }

            environmentalContexts.add(ImmutableEnvironmentalContext.builder()
                    .term(environmentContextObject.has("term") ? environmentContextObject.getAsJsonPrimitive("term").getAsString() : null)
                    .description(environmentContextObject.getAsJsonPrimitive("description").getAsString())
                    .taxonomy(environmentContextObject.has("taxonomy")
                            ? createTaxonomy(environmentContextObject.getAsJsonObject("taxonomy"))
                            : null)
                    .source(environmentContextObject.has("source")
                            ? environmentContextObject.getAsJsonPrimitive("source").getAsString()
                            : null)
                    .usanStem(environmentContextObject.has("usan_stem") ? environmentContextObject.getAsJsonPrimitive("usan_stem")
                            .getAsString() : null)
                    .approvedCountries(approvedCountries)
                    .id(environmentContextObject.has("id") && !environmentContextObject.get("id").isJsonNull()
                            ? environmentContextObject.getAsJsonPrimitive("id").getAsString()
                            : null)
                    .build());

        }
        return environmentalContexts;
    }

    @NotNull
    private static Taxonomy createTaxonomy(@NotNull JsonObject environmentContextObject) {
        return ImmutableTaxonomy.builder()
                .kingdom(environmentContextObject.getAsJsonPrimitive("kingdom").getAsString())
                .directParent(environmentContextObject.getAsJsonPrimitive("direct-parent").getAsString())
                .classs(environmentContextObject.getAsJsonPrimitive("class").getAsString())
                .subClass(environmentContextObject.has("subclass")
                        ? environmentContextObject.getAsJsonPrimitive("subclass").getAsString()
                        : null)
                .superClass(environmentContextObject.getAsJsonPrimitive("superclass").getAsString())
                .build();
    }

    @NotNull
    private static List<Evidence> createEvidence(@NotNull JsonArray evidenceArray) {
        List<Evidence> listEvidence = Lists.newArrayList();

        for (JsonElement evidenceElement : evidenceArray) {
            JsonObject evidenceObject = evidenceElement.getAsJsonObject();
            Set<String> keysEvidence = evidenceObject.keySet();
            if (!EXPECTED_EVIDENCE_ELEMENT_SIZES.contains(keysEvidence.size())) {
                LOGGER.warn("Found " + keysEvidence.size() + " in evidence rather than the expected " + EXPECTED_EVIDENCE_ELEMENT_SIZES);
                LOGGER.warn(keysEvidence);
            }

            listEvidence.add(ImmutableEvidence.builder()
                    .info(!evidenceObject.get("info").isJsonNull() ? createEvidenceInfo(evidenceObject.getAsJsonObject("info")) : null)
                    .evidenceType(createEvidenceType(evidenceObject.getAsJsonObject("evidenceType")))
                    .description(!evidenceObject.get("description").isJsonNull() ? evidenceObject.getAsJsonPrimitive("description")
                            .getAsString() : null)
                    .build());
        }
        return listEvidence;
    }

    @NotNull
    private static EvidenceType createEvidenceType(@NotNull JsonObject evidenceTypeObject) {
        if (!EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES.contains(evidenceTypeObject.keySet().size())) {
            LOGGER.warn("Found " + evidenceTypeObject.keySet().size() + " in evidence type rather than the expected "
                    + EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES);
            LOGGER.warn(evidenceTypeObject.keySet());
        }
        return ImmutableEvidenceType.builder()
                .sourceName(evidenceTypeObject.getAsJsonPrimitive("sourceName").getAsString())
                .id(evidenceTypeObject.has("id") ? evidenceTypeObject.getAsJsonPrimitive("id").getAsString() : null)
                .build();
    }

    @NotNull
    private static EvidenceInfo createEvidenceInfo(@NotNull JsonObject evidenceInfoObject) {
        if (!EXPECTED_EVIDENCE_INFO_ELEMENT_SIZES.contains(evidenceInfoObject.keySet().size())) {
            LOGGER.warn("Found " + evidenceInfoObject.keySet().size() + " in evidence info rather than the expected "
                    + EXPECTED_EVIDENCE_INFO_ELEMENT_SIZES);
            LOGGER.warn(evidenceInfoObject.keySet());
        }
        return ImmutableEvidenceInfo.builder()
                .publications(jsonArrayToStringList(evidenceInfoObject.getAsJsonArray("publications")))
                .build();
    }

    @NotNull
    private static Phenotype createPhenotype(@NotNull JsonObject phenotypeObject) {
        if (!EXPECTED_PHENOTYPE_ELEMENT_SIZES.contains(phenotypeObject.keySet().size())) {
            LOGGER.warn("Found " + phenotypeObject.keySet().size() + " in phenotype rather than the expected "
                    + EXPECTED_PHENOTYPE_ELEMENT_SIZES);
            LOGGER.warn(phenotypeObject.keySet());
        }
        return ImmutablePhenotype.builder()
                .type(phenotypeObject.has("type") ? createPhenotypeType(phenotypeObject.getAsJsonObject("type")) : null)
                .description(phenotypeObject.getAsJsonPrimitive("description").getAsString())
                .family(phenotypeObject.getAsJsonPrimitive("family").getAsString())
                .id(phenotypeObject.has("id") ? phenotypeObject.getAsJsonPrimitive("id").getAsString() : null)
                .build();
    }

    @NotNull
    private static PhenotypeType createPhenotypeType(JsonObject phenotypeTypeObject) {
        if (!EXPECTED_PHENOTYPE_TYPE_ELEMENT_SIZES.contains(phenotypeTypeObject.keySet().size())) {
            LOGGER.warn("Found " + phenotypeTypeObject.keySet().size() + " in phenotype type rather than the expected "
                    + EXPECTED_PHENOTYPE_TYPE_ELEMENT_SIZES);
            LOGGER.warn(phenotypeTypeObject.keySet());
        }
        return ImmutablePhenotypeType.builder()
                .source(!phenotypeTypeObject.get("source").isJsonNull()
                        ? phenotypeTypeObject.getAsJsonPrimitive("source").getAsString()
                        : null)
                .term(phenotypeTypeObject.getAsJsonPrimitive("term").getAsString())
                .id(phenotypeTypeObject.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static List<String> jsonArrayToStringList(@NotNull JsonArray array) {
        List<String> values = Lists.newArrayList();
        for (JsonElement element : array) {
            values.add(element.getAsString());
        }
        return values;
    }
}
