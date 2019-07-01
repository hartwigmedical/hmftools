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
import com.google.gson.annotations.JsonAdapter;
import com.google.gson.stream.JsonReader;
import com.google.gson.stream.JsonToken;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.BRCA;
import com.hartwig.hmftools.vicc.datamodel.BRCApart1;
import com.hartwig.hmftools.vicc.datamodel.BRCApart2;
import com.hartwig.hmftools.vicc.datamodel.BiologicalOncoKb;
import com.hartwig.hmftools.vicc.datamodel.Cgi;
import com.hartwig.hmftools.vicc.datamodel.ClinicalOncoKb;
import com.hartwig.hmftools.vicc.datamodel.ConsequenceOncoKb;
import com.hartwig.hmftools.vicc.datamodel.DrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.GeneOncokb;
import com.hartwig.hmftools.vicc.datamodel.GenePmkb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCA;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCApart1;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBRCApart2;
import com.hartwig.hmftools.vicc.datamodel.ImmutableBiologicalOncoKb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableCgi;
import com.hartwig.hmftools.vicc.datamodel.ImmutableClinicalOncoKb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableConsequenceOncoKb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableDrugAbstracts;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneOncokb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGenePmkb;
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
import com.hartwig.hmftools.vicc.datamodel.ImmutableTissuePmkb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTumorPmkb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableVariantOncokb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableVariantPmkb;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Jax;
import com.hartwig.hmftools.vicc.datamodel.JaxIndications;
import com.hartwig.hmftools.vicc.datamodel.JaxMolecularProfile;
import com.hartwig.hmftools.vicc.datamodel.JaxReferences;
import com.hartwig.hmftools.vicc.datamodel.JaxTherapy;
import com.hartwig.hmftools.vicc.datamodel.KbSpecificObject;
import com.hartwig.hmftools.vicc.datamodel.Oncokb;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Pmkb;
import com.hartwig.hmftools.vicc.datamodel.Sage;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.TissuePmkb;
import com.hartwig.hmftools.vicc.datamodel.TumorPmkb;
import com.hartwig.hmftools.vicc.datamodel.VariantOncokb;
import com.hartwig.hmftools.vicc.datamodel.VariantPmkb;
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

    private static final List<Integer> EXPECTED_CLINICAL_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(11);
    private static final List<Integer> EXPECTED_DRUGS_ABSTRACT_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(2);
    private static final List<Integer> EXPECTED_BIOLOGICAL_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(9);
    private static final List<Integer> EXPECTED_VARIANT_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(8);
    private static final List<Integer> EXPECTED_CONSEQUENCE_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(3);
    private static final List<Integer> EXPECTED_GENE_ONCOKB_ELEMENT_SIZES = Lists.newArrayList(8);

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
            }
            entries.add(viccEntryBuilder.build());

        }
        reader.close();

        return entries;
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
        return ImmutableJaxMolecularProfile.builder()
                .profileName(objectMolecularProfile.getAsJsonPrimitive("profileName").getAsString())
                .id(objectMolecularProfile.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    @NotNull
    private static JaxTherapy createJaxTherapy(@NotNull JsonObject objectTherapy) {
        return ImmutableJaxTherapy.builder()
                .id(objectTherapy.getAsJsonPrimitive("id").getAsString())
                .therapyName(objectTherapy.getAsJsonPrimitive("therapyName").getAsString())
                .build();
    }

    @NotNull
    private static JaxIndications createJaxIndications(@NotNull JsonObject objectIndications) {
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
            listReferences.add(ImmutableJaxReferences.builder()
                    .url(references.getAsJsonObject().getAsJsonPrimitive("url").getAsString())
                    .id(references.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .pubMedId(references.getAsJsonObject().get("pubMedId").isJsonNull() ? null :  references.getAsJsonObject().getAsJsonPrimitive("pubMedId").getAsString())
                    .title(references.getAsJsonObject().getAsJsonPrimitive("title").getAsString())
                    .build());
        }
        return listReferences;
    }

    @NotNull
    private static Oncokb createOncoKbBiological(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().biologicalOncoKb(createBiologicalOncoKb(objectOncoKb.getAsJsonObject("biological"))).build();
    }

    @NotNull
    private static Oncokb createOncoKbClinical(@NotNull JsonObject objectOncoKb) {
        return ImmutableOncokb.builder().clinicalOncoKb(createClinicalOncoKb(objectOncoKb.getAsJsonObject("clinical"))).build();
    }

    @NotNull
    private static ClinicalOncoKb createClinicalOncoKb(@NotNull JsonObject objectClinical) {
        Set<String> keysClinical = objectClinical.keySet();
        if (!EXPECTED_CLINICAL_ONCOKB_ELEMENT_SIZES.contains(keysClinical.size())) {
            LOGGER.warn("Found " + keysClinical.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_CLINICAL_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysClinical);
        }
        return ImmutableClinicalOncoKb.builder()
                .RefSeq(objectClinical.getAsJsonPrimitive("RefSeq").getAsString())
                .level(objectClinical.getAsJsonPrimitive("level").getAsString())
                .Isoform(objectClinical.getAsJsonPrimitive("Isoform").getAsString())
                .variantOncokb(createVariantOncoKb(objectClinical.getAsJsonObject("variant")))
                .entrezGeneID(objectClinical.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .drugPmids(objectClinical.getAsJsonPrimitive("drugPmids").getAsString())
                .cancerType(objectClinical.getAsJsonPrimitive("cancerType").getAsString())
                .drug(objectClinical.getAsJsonPrimitive("drug").getAsString())
                .gene(objectClinical.getAsJsonPrimitive("gene").getAsString())
                .levelLabel(objectClinical.getAsJsonPrimitive("level_label").getAsString())
                .drugAbstracts(createDrugsAbstracts(objectClinical.getAsJsonArray("drugAbstracts")))
                .build();
    }

    @NotNull
    private static List<DrugAbstracts> createDrugsAbstracts(@NotNull JsonArray arrayDrugsAbstracts) {

        List<DrugAbstracts> listDrugsabstracts = Lists.newArrayList();
        for (JsonElement drugAbstracts : arrayDrugsAbstracts) {
            Set<String> keysBiological = drugAbstracts.getAsJsonObject().keySet();

            if (!EXPECTED_DRUGS_ABSTRACT_ONCOKB_ELEMENT_SIZES.contains(keysBiological.size())) {
                LOGGER.warn("Found " + keysBiological.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_DRUGS_ABSTRACT_ONCOKB_ELEMENT_SIZES);
                LOGGER.warn(keysBiological);
            }
            listDrugsabstracts.add(ImmutableDrugAbstracts.builder()
                    .text(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("text").getAsString())
                    .link(drugAbstracts.getAsJsonObject().getAsJsonPrimitive("link").getAsString())
                    .build());
        }
        return listDrugsabstracts;
    }

    @NotNull
    private static BiologicalOncoKb createBiologicalOncoKb(@NotNull JsonObject objectBiological) {
        Set<String> keysBiological = objectBiological.keySet();
        if (!EXPECTED_BIOLOGICAL_ONCOKB_ELEMENT_SIZES.contains(keysBiological.size())) {
            LOGGER.warn("Found " + keysBiological.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_BIOLOGICAL_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysBiological);
        }

        return ImmutableBiologicalOncoKb.builder()
                .mutationEffectPmids(objectBiological.getAsJsonPrimitive("mutationEffectPmids").getAsString())
                .Isoform(objectBiological.getAsJsonPrimitive("Isoform").getAsString())
                .variantOncokb(createVariantOncoKb(objectBiological.getAsJsonObject("variant")))
                .entrezGeneID(objectBiological.getAsJsonPrimitive("Entrez Gene ID").getAsString())
                .oncogenic(objectBiological.getAsJsonPrimitive("oncogenic").getAsString())
                .mutationEffect(objectBiological.getAsJsonPrimitive("mutationEffect").getAsString())
                .RefSeq(objectBiological.getAsJsonPrimitive("RefSeq").getAsString())
                .gene(objectBiological.getAsJsonPrimitive("gene").getAsString())
                .mutationEffectAbstracts(objectBiological.getAsJsonPrimitive("mutationEffectAbstracts").getAsString())
                .build();
    }

    @NotNull
    private static VariantOncokb createVariantOncoKb(@NotNull JsonObject objectVariant) {
        Set<String> keysVariant = objectVariant.keySet();

        if (!EXPECTED_VARIANT_ONCOKB_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn("Found " + keysVariant.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_VARIANT_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }

        return ImmutableVariantOncokb.builder()
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
                .consequenceOncoKb(createConsequenceOncokb(objectVariant.getAsJsonObject("consequence")))
                .geneOncokb(createGeneOncoKb(objectVariant.getAsJsonObject("gene")))
                .build();
    }

    @NotNull
    private static ConsequenceOncoKb createConsequenceOncokb(@NotNull JsonObject objectConsequence) {
        Set<String> keysConsequence = objectConsequence.keySet();

        if (!EXPECTED_CONSEQUENCE_ONCOKB_ELEMENT_SIZES.contains(keysConsequence.size())) {
            LOGGER.warn("Found " + keysConsequence.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_CONSEQUENCE_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysConsequence);
        }

        return ImmutableConsequenceOncoKb.builder()
                .term(objectConsequence.getAsJsonPrimitive("term").getAsString())
                .description(objectConsequence.getAsJsonPrimitive("description").getAsString())
                .isGenerallyTruncating(objectConsequence.getAsJsonPrimitive("isGenerallyTruncating").getAsString())
                .build();
    }

    @NotNull
    private static GeneOncokb createGeneOncoKb(@NotNull JsonObject objectGene) {
        Set<String> keysGene = objectGene.keySet();

        if (!EXPECTED_GENE_ONCOKB_ELEMENT_SIZES.contains(keysGene.size())) {
            LOGGER.warn("Found " + keysGene.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_GENE_ONCOKB_ELEMENT_SIZES);
            LOGGER.warn(keysGene);
        }

        return ImmutableGeneOncokb.builder()
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
    private static List<TumorPmkb> createTumor(@NotNull JsonObject tumor) {
        Set<String> keysTumor = tumor.keySet();
        if (!EXPECTED_PMKB_TUMOR_ELEMENT_SIZES.contains(keysTumor.size())) {
            LOGGER.warn("Found " + keysTumor.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_PMKB_TUMOR_ELEMENT_SIZES);
            LOGGER.warn(keysTumor);
        }

        List<TumorPmkb> listTumor = Lists.newArrayList();
        listTumor.add(ImmutableTumorPmkb.builder()
                .id(tumor.getAsJsonPrimitive("id").getAsString())
                .name(tumor.getAsJsonPrimitive("name").getAsString())
                .build());

        return listTumor;
    }

    @NotNull
    private static List<TissuePmkb> createTissue(@NotNull JsonArray tissues) {

        List<TissuePmkb> listTissue = Lists.newArrayList();
        for (JsonElement tissue : tissues) {
            Set<String> keysTissue = tissue.getAsJsonObject().keySet();
            if (!EXPECTED_PMKB_TISSUE_ELEMENT_SIZES.contains(keysTissue.size())) {
                LOGGER.warn("Found " + keysTissue.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_PMKB_TISSUE_ELEMENT_SIZES);
                LOGGER.warn(keysTissue);
            }
            listTissue.add(ImmutableTissuePmkb.builder()
                    .id(tissue.getAsJsonObject().getAsJsonPrimitive("id").getAsString())
                    .name(tissue.getAsJsonObject().getAsJsonPrimitive("name").getAsString())
                    .build());
        }
        return listTissue;
    }

    @NotNull
    private static List<VariantPmkb> createVariantPmkb(@NotNull JsonObject variant) {
        Set<String> keysVariant = variant.keySet();
        if (!EXPECTED_PMKB_VARIANT_ELEMENT_SIZES.contains(keysVariant.size())) {
            LOGGER.warn("Found " + keysVariant.size() + " elements in a vicc entry rather than the expected "
                    + EXPECTED_PMKB_VARIANT_ELEMENT_SIZES);
            LOGGER.warn(keysVariant);
        }
        return Lists.newArrayList(ImmutableVariantPmkb.builder()
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
    private static List<GenePmkb> createGene(@NotNull JsonObject variant) {
        JsonObject gene = variant.getAsJsonObject("gene");

        Set<String> keysgene = gene.keySet();
        if (!EXPECTED_PMKB_GENE_ELEMENT_SIZES.contains(keysgene.size())) {
            LOGGER.warn(
                    "Found " + keysgene.size() + " elements in a vicc entry rather than the expected " + EXPECTED_PMKB_GENE_ELEMENT_SIZES);
            LOGGER.warn(keysgene);
        }

        List<GenePmkb> listGene = Lists.newArrayList();
        listGene.add(ImmutableGenePmkb.builder()
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
