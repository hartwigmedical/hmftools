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

    // SAGE records hold 8 field (no "feature names") while all other knowledgebases hold 9 records.
    private static final List<Integer> EXPECTED_VICC_ENTRY_SIZES = Lists.newArrayList(8, 9);

    private static final List<Integer> EXPECTED_ASSOCIATION_ELEMENT_SIZES = Lists.newArrayList(9, 10);

    private ViccJsonReader() {
    }

    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        List<ViccEntry> entries = Lists.newArrayList();
        LOGGER.info("Reading VICC knowledgebase from " + jsonPath);
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject viccEntryElement = parser.parse(reader).getAsJsonObject();

            if (!EXPECTED_VICC_ENTRY_SIZES.contains(viccEntryElement.size())) {
                LOGGER.warn("Found " + viccEntryElement.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_VICC_ENTRY_SIZES);
                LOGGER.warn(viccEntryElement);
            }

            ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
            viccEntryBuilder.source(viccEntryElement.getAsJsonPrimitive("source").getAsString());
            viccEntryBuilder.genes(jsonArrayToStringList(viccEntryElement.getAsJsonArray("genes")));

            viccEntryBuilder.geneIdentifiers(createGeneIdentifiers(viccEntryElement));

            if (viccEntryElement.has("feature_names")) {
                JsonElement featureNames = viccEntryElement.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(createFeatures());

            JsonElement elementAssociation = viccEntryElement.get("association");
            Set<String> keysAssociation = elementAssociation.getAsJsonObject().keySet();

//            if (!EXPECTED_ASSOCIATION_ELEMENT_SIZES.contains(keysAssociation.size())) {
//                LOGGER.warn("Found " + keysAssociation.size() + " elements in a vicc entry rather than the expected "
//                        + EXPECTED_ASSOCIATION_ELEMENT_SIZES);
//                LOGGER.warn(keysAssociation);
//            }

            viccEntryBuilder.association(createAssociationEmpty());

            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("dev_tags")));

            entries.add(viccEntryBuilder.build());
        }
        reader.close();

        return entries;
    }

    @NotNull
    private static List<Feature> createFeatures() {
        List<Feature> featureList = Lists.newArrayList();
        featureList.add(ImmutableFeature.builder()
                .name(Strings.EMPTY)
                .biomarkerType(Strings.EMPTY)
                .referenceName(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .start(Strings.EMPTY)
                .end(Strings.EMPTY)
                .ref(Strings.EMPTY)
                .alt(Strings.EMPTY)
                .provenance(Lists.newArrayList())
                .provenanceRule(Strings.EMPTY)
                .geneSymbol(Strings.EMPTY)
                .synonyms(Lists.newArrayList())
                .entrezId(Strings.EMPTY)
                .sequenceOntology(createSequenceOntology())
                .links(Lists.newArrayList())
                .description(Strings.EMPTY)
                .build());
        return featureList;
    }

    @NotNull
    private static SequenceOntology createSequenceOntology() {
        return ImmutableSequenceOntology.builder()
                .hierarchy(Lists.newArrayList())
                .soid(Strings.EMPTY)
                .parentSoid(Strings.EMPTY)
                .name(Strings.EMPTY)
                .parentName(Strings.EMPTY)
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

    @NotNull
    private static List<GeneIdentifier> createGeneIdentifiers(JsonObject viccEntryElement) {
        JsonArray geneIdentifiers = viccEntryElement.getAsJsonArray("gene_identifiers");
        List<GeneIdentifier> listGeneIdentifiers = Lists.newArrayList();

        for (JsonElement elementGeneIdentifier : geneIdentifiers) {
            listGeneIdentifiers.add(toGeneIdentifier(elementGeneIdentifier.getAsJsonObject()));

        }
        return listGeneIdentifiers;
    }

    @NotNull
    private static GeneIdentifier toGeneIdentifier(@NotNull JsonObject objectGeneIdentifier) {
        return ImmutableGeneIdentifier.builder()
                .symbol(objectGeneIdentifier.getAsJsonPrimitive("symbol").getAsString())
                .entrezId(objectGeneIdentifier.getAsJsonPrimitive("entrez_id").getAsString())
                .ensemblGeneId(Strings.EMPTY)
                .ensemblGeneId(Strings.EMPTY)
                .build();
    }

    @NotNull
    private static Association createAssociationEmpty() {
        return ImmutableAssociation.builder()
                .variantName(Strings.EMPTY)
                .evidence(Lists.newArrayList())
                .evidenceLevel(Strings.EMPTY)
                .evidenceLabel(Strings.EMPTY)
                .responseType(Strings.EMPTY)
                .drugLabels(Strings.EMPTY)
                .sourceLink(Strings.EMPTY)
                .publicationUrls(Lists.newArrayList())
                .phenotype(ImmutablePhenotype.builder()
                                .type(ImmutablePhenotypeType.builder()
                                        .source(Strings.EMPTY)
                                        .term(Strings.EMPTY)
                                        .id(Strings.EMPTY)
                                        .build())
                                .description(Strings.EMPTY)
                                .family(Strings.EMPTY)
                                .build())
                .description(Strings.EMPTY)
                .environmentalContexts(Lists.newArrayList())
                .oncogenic(Strings.EMPTY)
                .build();
    }

    @NotNull
    private static Association createAssociation(JsonObject objectAssociation) {
        return ImmutableAssociation.builder()
                .variantName(objectAssociation.has("variant_name")
                        ? objectAssociation.getAsJsonPrimitive("variant_name").getAsString()
                        : null)
                .evidence(createEvidence(objectAssociation.getAsJsonArray("evidence")))
                .evidenceLevel(objectAssociation.getAsJsonPrimitive("evidence_level").getAsString())
                .evidenceLabel(objectAssociation.getAsJsonPrimitive("evidence_label").getAsString())
                .responseType(objectAssociation.getAsJsonPrimitive("response_type").getAsString())
                .drugLabels(objectAssociation.getAsJsonPrimitive("drug_labels").getAsString())
                .sourceLink(objectAssociation.has("source_link") ? objectAssociation.getAsJsonPrimitive("source_link").getAsString() : null)
                .publicationUrls(Lists.newArrayList(objectAssociation.getAsJsonPrimitive("publication_url").getAsString()))
                .phenotype(createPhenotype(objectAssociation.getAsJsonObject("phenotype")))
                .description(objectAssociation.getAsJsonPrimitive("description").getAsString())
                .environmentalContexts(createEnvironmentalContexts(objectAssociation.getAsJsonArray("environmentalContexts")))
                .oncogenic(objectAssociation.has("oncogenic") ? objectAssociation.getAsJsonPrimitive("oncogenic").getAsString() : null)
                .build();
    }

    @NotNull
    private static List<EnvironmentalContext> createEnvironmentalContexts(@NotNull JsonArray arrayEnvironmentalContexts) {
        List<EnvironmentalContext> environmentalContexts = Lists.newArrayList();

        for (JsonElement elementEnvironmentContext : arrayEnvironmentalContexts) {
            JsonObject objectEnvironmentContext = elementEnvironmentContext.getAsJsonObject();

            List<String> approvedCountries = Lists.newArrayList();
            if (objectEnvironmentContext.has("approved_countries")) {
                for (JsonElement approvedCountriesElement : objectEnvironmentContext.getAsJsonArray("approved_countries")) {
                    approvedCountries.add(approvedCountriesElement.getAsString());
                }
            }

            environmentalContexts.add(ImmutableEnvironmentalContext.builder()
                    .term(objectEnvironmentContext.has("term") ? objectEnvironmentContext.getAsJsonPrimitive("term").getAsString() : null)
                    .description(objectEnvironmentContext.getAsJsonPrimitive("description").getAsString())
                    .taxonomy(objectEnvironmentContext.has("taxonomy")
                            ? createTaxonomy(objectEnvironmentContext.getAsJsonObject("taxonomy"))
                            : null)
                    .source(objectEnvironmentContext.has("source")
                            ? objectEnvironmentContext.getAsJsonPrimitive("source").getAsString()
                            : null)
                    .usanStem(objectEnvironmentContext.has("usan_stem") ? objectEnvironmentContext.getAsJsonPrimitive("usan_stem")
                            .getAsString() : null)
                    .approvedCountries(approvedCountries)
                    .id(objectEnvironmentContext.has("id") ? objectEnvironmentContext.getAsJsonPrimitive("id").getAsString() : null)
                    .build());

        }
        return environmentalContexts;
    }

    private static Taxonomy createTaxonomy(JsonObject elementEnvironmentContext) {
        return ImmutableTaxonomy.builder()
                .kingdom(elementEnvironmentContext.getAsJsonPrimitive("kingdom").getAsString())
                .directParent(elementEnvironmentContext.getAsJsonPrimitive("direct-parent").getAsString())
                .classs(elementEnvironmentContext.getAsJsonPrimitive("class").getAsString())
                .subClass(elementEnvironmentContext.has("subclass")
                        ? elementEnvironmentContext.getAsJsonPrimitive("subclass").getAsString()
                        : null)
                .superClass(elementEnvironmentContext.getAsJsonPrimitive("superclass").getAsString())
                .build();
    }

    private static List<Evidence> createEvidence(JsonArray arrayEvidence) {
        List<Evidence> listEvidence = Lists.newArrayList();

        for (JsonElement elementEvidence : arrayEvidence) {
            listEvidence.add(ImmutableEvidence.builder()
                    .info(createEvidenceInfo(elementEvidence.getAsJsonObject().getAsJsonObject("info")))
                    .evidenceType(createEvidenceType(elementEvidence.getAsJsonObject().getAsJsonObject("evidenceType")))
                    .description(elementEvidence.getAsJsonObject().getAsJsonPrimitive("description").getAsString())
                    .build());
        }
        return listEvidence;
    }

    private static EvidenceType createEvidenceType(JsonObject objectEvidenceType) {
        return ImmutableEvidenceType.builder()
                .sourceName(objectEvidenceType.getAsJsonPrimitive("sourceName").getAsString())
                .id(objectEvidenceType.has("id") ? objectEvidenceType.getAsJsonPrimitive("id").getAsString() : null)
                .build();
    }

    private static EvidenceInfo createEvidenceInfo(JsonObject objectEvidenceInfo) {
        List<String> publications = Lists.newArrayList();
        for (JsonElement publicationElements : objectEvidenceInfo.getAsJsonArray("publications")) {
            publications.add(publicationElements.getAsString());
        }

        return ImmutableEvidenceInfo.builder().publications(Lists.newArrayList(publications)).build();
    }

    private static PhenotypeType createPhenotypeType(JsonObject objectphenotypetype) {
        return ImmutablePhenotypeType.builder()
                .source(objectphenotypetype.getAsJsonPrimitive("source").getAsString())
                .term(objectphenotypetype.getAsJsonPrimitive("term").getAsString())
                .id(objectphenotypetype.getAsJsonPrimitive("id").getAsString())
                .build();
    }

    private static Phenotype createPhenotype(JsonObject objectPehnotype) {
        return ImmutablePhenotype.builder()
                .type(createPhenotypeType(objectPehnotype.getAsJsonObject("type")))
                .description(objectPehnotype.getAsJsonPrimitive("description").getAsString())
                .family(objectPehnotype.getAsJsonPrimitive("family").getAsString())
                .build();
    }
}
