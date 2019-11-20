package com.hartwig.hmftools.vicc.reader;

import static com.hartwig.hmftools.vicc.reader.JsonFunctions.jsonArrayToStringList;

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
import com.hartwig.hmftools.vicc.datamodel.EnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.Evidence;
import com.hartwig.hmftools.vicc.datamodel.EvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.EvidenceType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableSequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTaxonomy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.SequenceOntology;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccJsonReader {

    private static final Logger LOGGER = LogManager.getLogger(ViccJsonReader.class);

    private static final List<Integer> EXPECTED_EVIDENCE_INFO_ELEMENT_SIZES = Lists.newArrayList(1);
    private static final List<Integer> EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES = Lists.newArrayList(1, 2);
    private static final List<Integer> EXPECTED_PHENOTYPE_ELEMENT_SIZES = Lists.newArrayList(2, 3, 4);
    private static final List<Integer> EXPECTED_PHENOTYPE_TYPE_ELEMENT_SIZES = Lists.newArrayList(3);

    private static final List<Integer> EXPECTED_TAXONOMY = Lists.newArrayList(4, 5);
    private static final List<Integer> EXPECTED_ENVIRONMENT_CONTEXT = Lists.newArrayList(2, 4, 5, 6, 7, 8);

    private ViccJsonReader() {
    }

    @NotNull
    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        List<ViccEntry> entries = Lists.newArrayList();
        LOGGER.info("Reading VICC knowledgebase from {}", jsonPath);

        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryObject = parser.parse(reader).getAsJsonObject();
            ViccDatamodelCheckerFactory.viccEntryDatamodelChecker().check(viccEntryObject);

            ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
            viccEntryBuilder.source(viccEntryObject.getAsJsonPrimitive("source").getAsString());
            viccEntryBuilder.genes(jsonArrayToStringList(viccEntryObject.getAsJsonArray("genes")));
            viccEntryBuilder.geneIdentifiers(createGeneIdentifiers(viccEntryObject.getAsJsonArray("gene_identifiers")));

            // SAGE records have no "feature names" while all other knowledgebases do have it.
            if (viccEntryObject.has("feature_names")) {
                JsonElement featureNames = viccEntryObject.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(createFeatures(viccEntryObject.getAsJsonArray("features")));
            viccEntryBuilder.association(createAssociation(viccEntryObject.getAsJsonObject("association")));
            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryObject.getAsJsonArray("dev_tags")));

            if (viccEntryObject.has("cgi")) {
                viccEntryBuilder.KbSpecificObject(CgiObjectFactory.create(viccEntryObject.getAsJsonObject("cgi")));
            } else if (viccEntryObject.has("brca")) {
                viccEntryBuilder.KbSpecificObject(BRCAObjectFactory.create(viccEntryObject.getAsJsonObject("brca")));
            } else if (viccEntryObject.has("sage")) {
                viccEntryBuilder.KbSpecificObject(SageObjectFactory.create(viccEntryObject.getAsJsonObject("sage")));
            } else if (viccEntryObject.has("pmkb")) {
                viccEntryBuilder.KbSpecificObject(PmkbObjectFactory.create(viccEntryObject.getAsJsonObject("pmkb")));
            } else if (viccEntryObject.has("oncokb")) {
                viccEntryBuilder.KbSpecificObject(OncokbObjectFactory.create(viccEntryObject.getAsJsonObject("oncokb")));
            } else if (viccEntryObject.has("jax")) {
                viccEntryBuilder.KbSpecificObject(JaxObjectFactory.create(viccEntryObject.getAsJsonObject("jax")));
            } else if (viccEntryObject.has("jax_trials")) {
                viccEntryBuilder.KbSpecificObject(JaxTrialsObjectFactory.create(viccEntryObject.getAsJsonObject("jax_trials")));
            } else if (viccEntryObject.has("molecularmatch")) {
                viccEntryBuilder.KbSpecificObject(MolecularMatchObjectFactory.create(viccEntryObject.getAsJsonObject("molecularmatch")));
            } else if (viccEntryObject.has("molecularmatch_trials")) {
                viccEntryBuilder.KbSpecificObject(MolecularMatchTrialsObjectFactory.create(viccEntryObject.getAsJsonObject(
                        "molecularmatch_trials")));
            } else if (viccEntryObject.has("civic")) {
                viccEntryBuilder.KbSpecificObject(CivicObjectFactory.create(viccEntryObject.getAsJsonObject("civic")));
            } else {
                LOGGER.warn("Could not resolve kb specific object for {}", viccEntryObject);
            }

            entries.add(viccEntryBuilder.build());
        }

        reader.close();

        return entries;
    }

    @NotNull
    private static List<GeneIdentifier> createGeneIdentifiers(@NotNull JsonArray geneIdentifierArray) {
        List<GeneIdentifier> geneIdentifierList = Lists.newArrayList();

        for (JsonElement geneIdentifierElement : geneIdentifierArray) {
            JsonObject geneIdentifierObject = geneIdentifierElement.getAsJsonObject();
            ViccDatamodelCheckerFactory.geneIdentifierDatamodelChecker().check(geneIdentifierObject);

            geneIdentifierList.add(ImmutableGeneIdentifier.builder()
                    .symbol(geneIdentifierObject.getAsJsonPrimitive("symbol").getAsString())
                    .entrezId(geneIdentifierObject.getAsJsonPrimitive("entrez_id").getAsString())
                    .ensemblGeneId(!geneIdentifierObject.get("ensembl_gene_id").isJsonNull() ? geneIdentifierObject.getAsJsonPrimitive(
                            "ensembl_gene_id").getAsString() : null)
                    .build());
        }

        return geneIdentifierList;
    }

    @NotNull
    private static List<Feature> createFeatures(@NotNull JsonArray featureArray) {
        List<Feature> featureList = Lists.newArrayList();
        ViccDatamodelChecker featureDatamodelChecker = ViccDatamodelCheckerFactory.featureDatamodelChecker();

        for (JsonElement featureElement : featureArray) {
            JsonObject featureObject = featureElement.getAsJsonObject();
            featureDatamodelChecker.check(featureObject);

            featureList.add(ImmutableFeature.builder()
                    .name(featureObject.has("name") ? featureObject.getAsJsonPrimitive("name").getAsString() : null)
                    .biomarkerType(featureObject.has("biomarker_type")
                            ? featureObject.getAsJsonPrimitive("biomarker_type").getAsString()
                            : null)
                    .referenceName(featureObject.has("referenceName")
                            ? featureObject.getAsJsonPrimitive("referenceName").getAsString()
                            : null)
                    .chromosome(featureObject.has("chromosome") ? featureObject.getAsJsonPrimitive("chromosome").getAsString() : null)
                    .start(featureObject.has("start") && !featureObject.get("start").isJsonNull()
                            ? featureObject.getAsJsonPrimitive("start").getAsString()
                            : null)
                    .end(featureObject.has("end") && !featureObject.get("end").isJsonNull() ? featureObject.getAsJsonPrimitive("end")
                            .getAsString() : null)
                    .ref(featureObject.has("ref") && !featureObject.get("ref").isJsonNull() ? featureObject.getAsJsonPrimitive("ref")
                            .getAsString() : null)
                    .alt(featureObject.has("alt") && !featureObject.get("alt").isJsonNull() ? featureObject.getAsJsonPrimitive("alt")
                            .getAsString() : null)
                    .provenance(Lists.newArrayList())
                    .provenanceRule(featureObject.has("provenance_rule")
                            ? featureObject.getAsJsonPrimitive("provenance_rule").getAsString()
                            : null)
                    .geneSymbol(featureObject.has("geneSymbol") && !featureObject.get("geneSymbol").isJsonNull()
                            ? featureObject.getAsJsonPrimitive("geneSymbol").getAsString()
                            : null)
                    .synonyms(featureObject.has("synonyms") ? jsonArrayToStringList(featureObject.getAsJsonArray("synonyms")) : null)
                    .entrezId(featureObject.has("entrez_id") ? featureObject.getAsJsonPrimitive("entrez_id").getAsString() : null)
                    .sequenceOntology(featureObject.has("sequence_ontology") ? createSequenceOntology(featureObject.getAsJsonObject(
                            "sequence_ontology")) : null)
                    .links(featureObject.has("links") ? jsonArrayToStringList(featureObject.getAsJsonArray("links")) : null)
                    .description(featureObject.has("description") ? featureObject.getAsJsonPrimitive("description").getAsString() : null)
                    .build());
        }

        return featureList;
    }

    @NotNull
    private static SequenceOntology createSequenceOntology(@NotNull JsonObject objectSequenceOntology) {
        ViccDatamodelCheckerFactory.sequenceOntologyDatamodelChecker().check(objectSequenceOntology);

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
    private static Association createAssociation(@NotNull JsonObject associationObject) {
        ViccDatamodelCheckerFactory.associationDatamodelChecker().check(associationObject);

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
    private static List<Evidence> createEvidence(@NotNull JsonArray evidenceArray) {
        List<Evidence> evidenceList = Lists.newArrayList();
        ViccDatamodelChecker evidenceChecker = ViccDatamodelCheckerFactory.evidenceDatamodelChecker();

        for (JsonElement evidenceElement : evidenceArray) {
            JsonObject evidenceObject = evidenceElement.getAsJsonObject();
            evidenceChecker.check(evidenceObject);
            evidenceList.add(ImmutableEvidence.builder()
                    .info(!evidenceObject.get("info").isJsonNull() ? createEvidenceInfo(evidenceObject.getAsJsonObject("info")) : null)
                    .evidenceType(createEvidenceType(evidenceObject.getAsJsonObject("evidenceType")))
                    .description(!evidenceObject.get("description").isJsonNull() ? evidenceObject.getAsJsonPrimitive("description")
                            .getAsString() : null)
                    .build());
        }
        return evidenceList;
    }

    @NotNull
    private static EvidenceType createEvidenceType(@NotNull JsonObject evidenceTypeObject) {
        if (!EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES.contains(evidenceTypeObject.keySet().size())) {
            LOGGER.warn("Found {} in evidence type rather than the expected {}",
                    evidenceTypeObject.keySet().size(),
                    EXPECTED_EVIDENCE_TYPE_ELEMENT_SIZES);
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
            LOGGER.warn("Found {} in evidence info rather than the expected {}",
                    evidenceInfoObject.keySet().size(),
                    EXPECTED_EVIDENCE_INFO_ELEMENT_SIZES);
            LOGGER.warn(evidenceInfoObject.keySet());
        }
        return ImmutableEvidenceInfo.builder()
                .publications(jsonArrayToStringList(evidenceInfoObject.getAsJsonArray("publications")))
                .build();
    }

    @NotNull
    private static List<EnvironmentalContext> createEnvironmentalContexts(@NotNull JsonArray arrayEnvironmentalContexts) {
        List<EnvironmentalContext> environmentalContexts = Lists.newArrayList();

        for (JsonElement elementEnvironmentContext : arrayEnvironmentalContexts) {
            JsonObject environmentContextObject = elementEnvironmentContext.getAsJsonObject();
            Set<String> keysEnvironmentContext = environmentContextObject.keySet();

            if (!EXPECTED_ENVIRONMENT_CONTEXT.contains(keysEnvironmentContext.size())) {
                LOGGER.warn("Found {} in environmental context rather than the expected {}",
                        keysEnvironmentContext.size(),
                        EXPECTED_ENVIRONMENT_CONTEXT);
                LOGGER.warn(keysEnvironmentContext);
            }

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
                    .toxicity(environmentContextObject.has("toxicity") ? environmentContextObject.getAsJsonPrimitive("toxicity")
                            .getAsString() : null)
                    .id(environmentContextObject.has("id") && !environmentContextObject.get("id").isJsonNull()
                            ? environmentContextObject.getAsJsonPrimitive("id").getAsString()
                            : null)
                    .build());
        }
        return environmentalContexts;
    }

    @NotNull
    private static Taxonomy createTaxonomy(@NotNull JsonObject environmentContextObject) {
        if (!EXPECTED_TAXONOMY.contains(environmentContextObject.keySet().size())) {
            LOGGER.warn("Found {} in taxonomy rather than the expected {}", environmentContextObject.keySet().size(), EXPECTED_TAXONOMY);
            LOGGER.warn(environmentContextObject.keySet());
        }

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
    private static Phenotype createPhenotype(@NotNull JsonObject phenotypeObject) {
        if (!EXPECTED_PHENOTYPE_ELEMENT_SIZES.contains(phenotypeObject.keySet().size())) {
            LOGGER.warn("Found {} in phenotype rather than the expected {}",
                    phenotypeObject.keySet().size(),
                    EXPECTED_PHENOTYPE_ELEMENT_SIZES);
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
            LOGGER.warn("Found {} in phenotype type rather than the expected {}",
                    phenotypeTypeObject.keySet().size(),
                    EXPECTED_PHENOTYPE_TYPE_ELEMENT_SIZES);
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
}
