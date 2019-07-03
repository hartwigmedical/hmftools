package com.hartwig.hmftools.vicc;

import java.io.FileReader;
import java.io.IOException;
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
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrials;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsIndications;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsTherapies;
import com.hartwig.hmftools.vicc.datamodel.ImmutableJaxTrialsVariantRequirementDetails;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatch;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.ImmutableMolecularMatchVariantInfo;
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
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAst;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstLeft;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchAstRight;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchClassification;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchCriteriaUnmet;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchGRch37Location;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchLocations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchMutations;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchPrevalence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchSource;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTags;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTherapeuticContext;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTierExplanation;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequence;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchTranscriptConsequencesGRCH37;
import com.hartwig.hmftools.vicc.datamodel.MolecularMatchVariantInfo;
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
import org.omg.CORBA.PRIVATE_MEMBER;

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

    private static final List<Integer> EXPECTED_MOLECULARMATCH_ELEMENT_SIZES = Lists.newArrayList(36);

    private ViccJsonReader() {
    }

    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        List<ViccEntry> entries = Lists.newArrayList();
        LOGGER.info("Reading VICC knowledgebase from " + jsonPath);
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryObject = parser.parse(reader).getAsJsonObject();

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
                LOGGER.warn("Found " + keysAssociation.size() + " elements in a vicc entry rather than the expected "
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
                    LOGGER.warn(
                            "Found " + keysCgi.size() + " elements in a vicc entry rather than the expected " + EXPECTED_CGI_ELEMENT_SIZES);
                    LOGGER.warn(keysCgi);
                }
            }

            JsonObject objectBRCA = viccEntryObject.getAsJsonObject("brca");
            if (viccEntryObject.has("brca")) {
                Set<String> keysBRCA = objectBRCA.keySet();
                if (!EXPECTED_BRCA_ELEMENT_SIZES.contains(keysBRCA.size())) {
                    LOGGER.warn("Found " + keysBRCA.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_BRCA_ELEMENT_SIZES);
                    LOGGER.warn(keysBRCA);
                }
            }

            JsonObject objectSage = viccEntryObject.getAsJsonObject("sage");
            if (viccEntryObject.has("sage")) {
                Set<String> keysSage = objectSage.keySet();
                if (!EXPECTED_SAGE_ELEMENT_SIZES.contains(keysSage.size())) {
                    LOGGER.warn("Found " + keysSage.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_SAGE_ELEMENT_SIZES);
                    LOGGER.warn(keysSage);
                }
            }

            JsonObject objectPmkb = viccEntryObject.getAsJsonObject("pmkb");
            if (viccEntryObject.has("pmkb")) {
                Set<String> keysPmkb = objectPmkb.keySet();
                if (!EXPECTED_PMKB_ELEMENT_SIZES.contains(keysPmkb.size())) {
                    LOGGER.warn("Found " + keysPmkb.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_PMKB_ELEMENT_SIZES);
                    LOGGER.warn(keysPmkb);
                }
            }

            JsonObject objectOncokb = viccEntryObject.getAsJsonObject("oncokb");
            if (viccEntryObject.has("oncokb")) {
                Set<String> keysOncokb = objectOncokb.keySet();
                if (!EXPECTED_ONCOKB_ELEMENT_SIZES.contains(keysOncokb.size())) {
                    LOGGER.warn("Found " + keysOncokb.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_ONCOKB_ELEMENT_SIZES);
                    LOGGER.warn(keysOncokb);
                }
            }

            JsonObject objectJax = viccEntryObject.getAsJsonObject("jax");
            if (viccEntryObject.has("jax")) {
                Set<String> keysJax = objectJax.keySet();
                if (!EXPECTED_JAX_ELEMENT_SIZES.contains(keysJax.size())) {
                    LOGGER.warn(
                            "Found " + keysJax.size() + " elements in a vicc entry rather than the expected " + EXPECTED_JAX_ELEMENT_SIZES);
                    LOGGER.warn(keysJax);
                }
            }

            JsonObject objectJaxTrials = viccEntryObject.getAsJsonObject("jax_trials");

            if (viccEntryObject.has("jax_trials")) {
                Set<String> keysJaxTrials = objectJaxTrials.keySet();
                if (!EXPECTED_JAX_TRIALS_ELEMENT_SIZES.contains(keysJaxTrials.size())) {
                    LOGGER.warn("Found " + keysJaxTrials.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_JAX_TRIALS_ELEMENT_SIZES);
                    LOGGER.warn(keysJaxTrials);
                }
            }

            JsonObject objectMolecularMatch = viccEntryObject.getAsJsonObject("molecularmatch");
            if (viccEntryObject.has("molecularmatch")) {
                Set<String> keysMolecularMatch = objectMolecularMatch.keySet();
                if (!EXPECTED_MOLECULARMATCH_ELEMENT_SIZES.contains(keysMolecularMatch.size())) {
                    LOGGER.warn("Found " + keysMolecularMatch.size() + " elements in a vicc entry rather than the expected "
                            + EXPECTED_MOLECULARMATCH_ELEMENT_SIZES);
                    LOGGER.warn(keysMolecularMatch);
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
            }
            entries.add(viccEntryBuilder.build());

        }
        reader.close();

        return entries;
    }

    @NotNull
    private static MolecularMatch createMolecularMatch(@NotNull JsonObject objectMolecularMatch) {
        return ImmutableMolecularMatch.builder()
                .criteriaUnmet(createCriteriaUnmet(objectMolecularMatch.getAsJsonArray("criteriaUnmet")))
                .prevalences(createPrevalence(objectMolecularMatch.getAsJsonArray("prevalence")))
                .score(objectMolecularMatch.getAsJsonPrimitive("_score").getAsString())
                .autoGenerateNarrative(objectMolecularMatch.getAsJsonPrimitive("autoGenerateNarrative").getAsString())
                .mutations(createMutations(objectMolecularMatch.getAsJsonArray("mutations")))
                .sources(createSource(objectMolecularMatch.getAsJsonArray("sources")))
                .clinicalSignificance(objectMolecularMatch.getAsJsonPrimitive("clinicalSignificance").getAsString())
                .id(objectMolecularMatch.getAsJsonPrimitive("id").getAsString())
                .includeCondition0(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeCondition0")))
                .includeCondition1(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeCondition1")))
                .uniqueKey(objectMolecularMatch.getAsJsonPrimitive("uniqueKey").getAsString())
                .civic(objectMolecularMatch.getAsJsonPrimitive("civic").getAsString())
                .hashKey(objectMolecularMatch.getAsJsonPrimitive("hashKey").getAsString())
                .regulatoryBodyApproved(objectMolecularMatch.getAsJsonPrimitive("regulatoryBodyApproved").getAsString())
                .version(objectMolecularMatch.getAsJsonPrimitive("version").getAsString())
                .includeMutation1(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeMutation1")))
                .guidelineBody(objectMolecularMatch.getAsJsonPrimitive("guidelineBody").getAsString())
                .regulatoryBody(objectMolecularMatch.getAsJsonPrimitive("regulatoryBody").getAsString())
                .customer(objectMolecularMatch.getAsJsonPrimitive("customer").getAsString())
                .direction(objectMolecularMatch.getAsJsonPrimitive("direction").getAsString())
                .ampcap(objectMolecularMatch.getAsJsonPrimitive("ampcap").getAsString())
                .asts(createAst(objectMolecularMatch.getAsJsonObject("ast")))
                .variantInfo(createVariantInfo(objectMolecularMatch.getAsJsonArray("variantInfo")))
                .tier(objectMolecularMatch.getAsJsonPrimitive("tier").getAsString())
                .tierExplanation(createTierExplanation(objectMolecularMatch.getAsJsonArray("tierExplanation")))
                .mvld(objectMolecularMatch.getAsJsonPrimitive("mvld").getAsString())
                .tags(createTags(objectMolecularMatch.getAsJsonArray("tags")))
                .criteriaMet(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("criteriaMet")))
                .biomarkerClass(objectMolecularMatch.getAsJsonPrimitive("biomarkerClass").getAsString())
                .classification(createClassification(objectMolecularMatch.getAsJsonArray("classifications")))
                .includeDrug1(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeDrug1")))
                .therapeuticContext(createTherapeuticContext(objectMolecularMatch.getAsJsonArray("therapeuticContext")))
                .sixtier(objectMolecularMatch.getAsJsonPrimitive("sixtier").getAsString())
                .narrative(objectMolecularMatch.getAsJsonPrimitive("narrative").getAsString())
                .expression(objectMolecularMatch.getAsJsonPrimitive("expression").getAsString())
                .includeGene0(jsonArrayToStringList(objectMolecularMatch.getAsJsonArray("includeGene0")))
                .build();
    }

    @NotNull
    private static MolecularMatchTherapeuticContext createTherapeuticContext(@NotNull JsonArray arrayTherapeuticContext) {
        MolecularMatchTherapeuticContext therapeuticContexts = ImmutableMolecularMatchTherapeuticContext.builder().build();
        for (JsonElement therapeuticContext: arrayTherapeuticContext) {
            therapeuticContexts = ImmutableMolecularMatchTherapeuticContext.builder()
                    .facet(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .suppress(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .valid(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .name(therapeuticContext.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build();
        }
        return therapeuticContexts;
    }

    @NotNull
    private static MolecularMatchClassification createClassification(@NotNull JsonArray objectClassifications) {
        MolecularMatchClassification classifications = ImmutableMolecularMatchClassification.builder().build();
        for (JsonElement classification: objectClassifications) {
            classifications = ImmutableMolecularMatchClassification.builder()
                    .end(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("End")))
                    .classification(classification.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .classificationOverride(classification.getAsJsonObject().getAsJsonPrimitive("classificationOverride").getAsString())
                    .start(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Start")))
                    .chr(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Chr")))
                    .geneSymbol(classification.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("pathology")))
                    .ref(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Ref")))
                    .description(classification.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .priority(classification.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .NucleotideChange(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("NucleotideChange")))
                    .parents(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("parents")))
                    .drugsExperimentalCount(classification.getAsJsonObject().getAsJsonPrimitive("drugsExperimentalCount").getAsString())
                    .exon(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Exon")))
                    .drugsApprovedOffLabelCount(classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOffLabelCount").getAsString())
                    .exonicFunc(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("ExonicFunc")))
                    .popFreqMax(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("PopFreqMax")))
                    .copyNumberType(classification.getAsJsonObject().getAsJsonPrimitive("copyNumberType").getAsString())
                    .publicationCount(classification.getAsJsonObject().getAsJsonPrimitive("publicationCount").getAsString())
                    .transcript(classification.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .dbSNP(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("dbSNP")))
                    .alt(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("Alt")))
                    .name(classification.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .rootTerm(classification.getAsJsonObject().getAsJsonPrimitive("rootTerm").getAsString())
                    .sources(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("sources")))
                    .drugsApprovedOnLabelCount(classification.getAsJsonObject().getAsJsonPrimitive("drugsApprovedOnLabelCount").getAsString())
                    .trialCount(classification.getAsJsonObject().getAsJsonPrimitive("trialCount").getAsString())
                    .alias(classification.getAsJsonObject().getAsJsonPrimitive("alias").getAsString())
                    .COSMIC_ID(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("COSMIC_ID")))
                    .transcripts(jsonArrayToStringList(classification.getAsJsonObject().getAsJsonArray("transcripts")))
                    .build();

        }
        return classifications;

    }
    @NotNull
    private static List<MolecularMatchTags> createTags(@NotNull JsonArray arrayTags) {
        List<MolecularMatchTags> tagsList = Lists.newArrayList();
        for (JsonElement tags: arrayTags) {
            tagsList.add(ImmutableMolecularMatchTags.builder()
                    .priority(tags.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(tags.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(tags.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(tags.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(tags.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(tags.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(tags.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(tags.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(tags.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .build());
        }
        return tagsList;
    }
    @NotNull
    private static List<MolecularMatchTierExplanation> createTierExplanation(@NotNull JsonArray arrarTierExplanation) {
        List<MolecularMatchTierExplanation> tierExplanationList = Lists.newArrayList();
        for (JsonElement tierExplanation: arrarTierExplanation) {
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

        for (JsonElement variantInfo: arrayVariantInfo){
            variantInfoList.add(ImmutableMolecularMatchVariantInfo.builder()
                    .classification(variantInfo.getAsJsonObject().getAsJsonPrimitive("classification").getAsString())
                    .name(variantInfo.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .consequences(jsonArrayToStringList(variantInfo.getAsJsonObject().getAsJsonArray("consequences")))
                    .fusions(jsonArrayToStringList(variantInfo.getAsJsonObject().getAsJsonArray("fusions")))
                    .locations(createLocations(variantInfo.getAsJsonObject().getAsJsonArray("locations")))
                    .geneFusionPartner(variantInfo.getAsJsonObject().getAsJsonPrimitive("geneFusionPartner").getAsString())
                    .COSMIC_ID(variantInfo.getAsJsonObject().getAsJsonPrimitive("COSMIC_ID").getAsString())
                    .gene(variantInfo.getAsJsonObject().getAsJsonPrimitive("gene").getAsString())
                    .transcript(variantInfo.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .popFreqMax(variantInfo.getAsJsonObject().getAsJsonPrimitive("popFreqMax").getAsString())
            .build());
        }
        return variantInfoList;
    }

    @NotNull
    private static List<MolecularMatchLocations> createLocations(@NotNull JsonArray arrayLocations) {
        List<MolecularMatchLocations> locationsList = Lists.newArrayList();
        for (JsonElement locations: arrayLocations) {
            locationsList.add(ImmutableMolecularMatchLocations.builder()
                    .aminoAcidChange(locations.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .intronNumber(locations.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(locations.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString())
                    .stop(locations.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(locations.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(locations.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(locations.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .alt(locations.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .referenceGenome(locations.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
                    .ref(locations.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .cdna(locations.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .build());
        }
        return locationsList;
    }
    @NotNull
    private static MolecularMatchAst createAst(@NotNull JsonObject objectAst) {
        return ImmutableMolecularMatchAst.builder()
                .operator(objectAst.getAsJsonPrimitive("operator").getAsString())
                .right(createRight(objectAst.getAsJsonObject("right")))
                .type(objectAst.getAsJsonPrimitive("type").getAsString())
                .left(createLeft(objectAst.getAsJsonObject("left")))
                .build();
    }

    @NotNull
    private static MolecularMatchAstLeft createLeft(@NotNull JsonObject objectLeft) {
        return ImmutableMolecularMatchAstLeft.builder()
                .raw(objectLeft.getAsJsonPrimitive("raw").getAsString())
                .type(objectLeft.getAsJsonPrimitive("type").getAsString())
                .value(objectLeft.getAsJsonPrimitive("value").getAsString())
                .build();
    }

    @NotNull
    private static MolecularMatchAstRight createRight(@NotNull JsonObject objectRight) {
        return ImmutableMolecularMatchAstRight.builder()
                .raw(objectRight.getAsJsonPrimitive("raw").getAsString())
                .type(objectRight.getAsJsonPrimitive("type").getAsString())
                .value(objectRight.getAsJsonPrimitive("value").getAsString())
                .build();

    }
    @NotNull
    private static List<MolecularMatchSource> createSource(@NotNull JsonArray arraySources) {
        List<MolecularMatchSource> sourcesList = Lists.newArrayList();
        for (JsonElement source : arraySources) {
            sourcesList.add(ImmutableMolecularMatchSource.builder()
                    .name(source.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .suppress(source.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .pubId(source.getAsJsonObject().getAsJsonPrimitive("pubId").getAsString())
                    .subType(source.getAsJsonObject().getAsJsonPrimitive("subType").getAsString())
                    .valid(source.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .link(source.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .year(source.getAsJsonObject().getAsJsonPrimitive("year").getAsString())
                    .trialId(source.getAsJsonObject().getAsJsonPrimitive("trialId").getAsString())
                    .type(source.getAsJsonObject().getAsJsonPrimitive("type").getAsString())
                    .id(source.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .build());
        }
        return sourcesList;
    }

    @NotNull
    private static List<MolecularMatchCriteriaUnmet> createCriteriaUnmet(@NotNull JsonArray arrayCriteriaUnmet) {
        List<MolecularMatchCriteriaUnmet> criteriaUnmetList = Lists.newArrayList();

        for (JsonElement criteriaUnmet : arrayCriteriaUnmet) {
            criteriaUnmetList.add(ImmutableMolecularMatchCriteriaUnmet.builder()
                    .priority(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("priority").getAsString())
                    .compositeKey(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .suppress(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .filterType(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("filterType").getAsString())
                    .term(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("term").getAsString())
                    .primary(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("primary").getAsString())
                    .facet(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("facet").getAsString())
                    .valid(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("valid").getAsString())
                    .custom(criteriaUnmet.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .build());
        }
        return criteriaUnmetList;
    }

    @NotNull
    private static List<MolecularMatchPrevalence> createPrevalence(@NotNull JsonArray arrayPrevelance) {
        List<MolecularMatchPrevalence> prevalenceList = Lists.newArrayList();

        for (JsonElement prevelance : arrayPrevelance) {
            prevalenceList.add(ImmutableMolecularMatchPrevalence.builder()
                    .count(prevelance.getAsJsonObject().getAsJsonPrimitive("count").getAsString())
                    .percent(prevelance.getAsJsonObject().getAsJsonPrimitive("percent").getAsString())
                    .studyId(prevelance.getAsJsonObject().getAsJsonPrimitive("studyId").getAsString())
                    .samples(prevelance.getAsJsonObject().getAsJsonPrimitive("samples").getAsString())
                    .molecular(prevelance.getAsJsonObject().getAsJsonPrimitive("molecular").getAsString())
                    .condition(prevelance.getAsJsonObject().getAsJsonPrimitive("condition").getAsString())
                    .build());
        }
        return prevalenceList;
    }

    @NotNull
    private static List<MolecularMatchMutations> createMutations(@NotNull JsonArray arrayMutations) {
        List<MolecularMatchMutations> mutationList = Lists.newArrayList();

        for (JsonElement mutation : arrayMutations) {
            mutationList.add(ImmutableMolecularMatchMutations.builder()
                    .transcriptConsequence(createTranscriptConsequence(mutation.getAsJsonObject().getAsJsonArray("transcriptConsequence")))
                    .longestTranscript(mutation.getAsJsonObject().getAsJsonPrimitive("longestTranscript").getAsString())
                    .description(mutation.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .mutationType(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("mutation_type")))
                    .src(mutation.getAsJsonObject().getAsJsonPrimitive("_src").getAsString())
                    .sources(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("sources")))
                    .synonyms(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("synonyms")))
                    .parents(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("parents")))
                    .gRch37Location(createGRCH37Location(mutation.getAsJsonObject().getAsJsonArray("GRCh37_location")))
                    .uniprotTranscript(mutation.getAsJsonObject().getAsJsonPrimitive("uniprotTranscript").getAsString())
                    .geneSymbol(mutation.getAsJsonObject().getAsJsonPrimitive("geneSymbol").getAsString())
                    .pathology(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("pathology")))
                    .transcript(mutation.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .id(mutation.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .cDNA(jsonArrayToStringList(mutation.getAsJsonObject().getAsJsonArray("cdna")))
                    .name(mutation.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());

        }
        return mutationList;

    }

    @NotNull
    private static MolecularMatchGRch37Location createGRCH37Location(@NotNull JsonArray arrayLocation) {
        MolecularMatchGRch37Location mutation = ImmutableMolecularMatchGRch37Location.builder().build();
        for (JsonElement location: arrayLocation) {
            mutation = ImmutableMolecularMatchGRch37Location.builder()
                    .compositeKey(location.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .ref(location.getAsJsonObject().getAsJsonPrimitive("ref").getAsString())
                    .stop(location.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .start(location.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(location.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .alt(location.getAsJsonObject().getAsJsonPrimitive("alt").getAsString())
                    .validated(location.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcriptConsequences(createConsequencesGRCH37(location.getAsJsonObject().getAsJsonArray("transcript_consequences")))
                    .strand(location.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .build();
        }
        return mutation;

    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequencesGRCH37> createConsequencesGRCH37(@NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequencesGRCH37> transcriptConsequencesGRCH37List = Lists.newArrayList();
        for (JsonElement transcriptConsequences: arrayTranscriptConsequence) {
            transcriptConsequencesGRCH37List.add(ImmutableMolecularMatchTranscriptConsequencesGRCH37.builder()
                    .aminoAcidChange(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .txSites(jsonArrayToStringList(transcriptConsequences.getAsJsonObject().getAsJsonArray("txSites")))
                    .exonNumber(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString())
                    .intronNumber(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .transcript(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(transcriptConsequences.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
            .build());
        }
        return transcriptConsequencesGRCH37List;
    }

    @NotNull
    private static List<MolecularMatchTranscriptConsequence> createTranscriptConsequence(@NotNull JsonArray arrayTranscriptConsequence) {
        List<MolecularMatchTranscriptConsequence> transcriptConsequenceList = Lists.newArrayList();

        for (JsonElement transcriptConsequence: arrayTranscriptConsequence) {
            transcriptConsequenceList.add(ImmutableMolecularMatchTranscriptConsequence.builder()
                    .aminoAcidChange(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("amino_acid_change").getAsString())
                    .compositeKey(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("compositeKey").getAsString())
                    .intronNumber(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("intronNumber").getAsString())
                    .exonNumber(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("exonNumber").getAsString())
                    .suppress(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("suppress").getAsString())
                    .stop(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("stop").getAsString())
                    .custom(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("custom").getAsString())
                    .start(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("start").getAsString())
                    .chr(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("chr").getAsString())
                    .strand(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("strand").getAsString())
                    .validated(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("validated").getAsString())
                    .transcript(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("transcript").getAsString())
                    .cdna(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("cdna").getAsString())
                    .referenceGenome(transcriptConsequence.getAsJsonObject().getAsJsonPrimitive("referenceGenome").getAsString())
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
                LOGGER.warn("Found " + keysIndications.size() + " elements in a vicc entry rather than the expected "
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
                LOGGER.warn("Found " + keysRequirementDetails.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysMolecularProfile.size() + " elements in a vicc entry rather than the expected "
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
                LOGGER.warn("Found " + keysTherapies.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysMolecularProfile.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysTherapy.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_JAX_THERAPY_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + keysIndications.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_JAX_INDICATIONS_SIZES);
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
                LOGGER.warn("Found " + keysReferences.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysClinical.size() + " elements in a vicc entry rather than the expected "
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
                .oncoKbDrugAbstracts(createDrugsAbstracts(objectClinical.getAsJsonArray("oncoKbDrugAbstracts")))
                .build();
    }

    @NotNull
    private static List<OncoKbDrugAbstracts> createDrugsAbstracts(@NotNull JsonArray arrayDrugsAbstracts) {

        List<OncoKbDrugAbstracts> listDrugsabstracts = Lists.newArrayList();
        for (JsonElement drugAbstracts : arrayDrugsAbstracts) {
            Set<String> keysBiological = drugAbstracts.getAsJsonObject().keySet();

            if (!EXPECTED_ONCOKB_DRUGS_ABSTRACT_ELEMENT_SIZES.contains(keysBiological.size())) {
                LOGGER.warn("Found " + keysBiological.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysBiological.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysVariant.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_ONCOKB_VARIANT_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + keysConsequence.size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + keysGene.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_ONCOKB_GENE_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + keysTumor.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_PMKB_TUMOR_ELEMENT_SIZES);
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
                LOGGER.warn("Found " + keysTissue.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_PMKB_TISSUE_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + keysVariant.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_PMKB_VARIANT_ELEMENT_SIZES);
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
            LOGGER.warn(
                    "Found " + keysgene.size() + " elements in a vicc entry rather than the expected " + EXPECTED_PMKB_GENE_ELEMENT_SIZES);
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
                LOGGER.warn("Found " + keysFeatures.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_FEATURES_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + keysSequenceOntology.size() + " elements in a vicc entry rather than the expected "
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
                LOGGER.warn("Found " + keysGeneIdentifier.size() + " elements in a vicc entry rather than the expected "
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
                LOGGER.warn("Found " + keysEvidence.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_EVIDENCE_ELEMENT_SIZES);
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
            LOGGER.warn("Found " + evidenceTypeObject.keySet().size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + evidenceInfoObject.keySet().size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + phenotypeObject.keySet().size() + " elements in a vicc entry rather than the expected "
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
            LOGGER.warn("Found " + phenotypeTypeObject.keySet().size() + " elements in a vicc entry rather than the expected "
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
