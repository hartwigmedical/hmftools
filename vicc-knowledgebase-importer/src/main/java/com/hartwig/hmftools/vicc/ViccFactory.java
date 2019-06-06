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
import com.hartwig.hmftools.vicc.datamodel.GeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEnvironmentalContext;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidence;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceInfo;
import com.hartwig.hmftools.vicc.datamodel.ImmutableEvidenceType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableGeneIdentifier;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableTaxonomy;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
import com.hartwig.hmftools.vicc.datamodel.PhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.Taxonomy;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class ViccFactory {
    private static final Logger LOGGER = LogManager.getLogger(ViccFactory.class);

    // SAGE records hold 8 field (no "feature names") while all other knowledgebases hold 9 records.
    private static final List<Integer> EXPECTED_VICC_ENTRY_SIZES = Lists.newArrayList(8, 9);

    private static final List<Integer> EXPECTED_ASSOCIATION_ELEMENT_SIZES = Lists.newArrayList(9, 10);

    private ViccFactory() {
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

            viccEntryBuilder.features(Lists.newArrayList());

            JsonElement elementAssociation = viccEntryElement.get("association");
            Set<String> keysAssociation = elementAssociation.getAsJsonObject().keySet();

            if (!EXPECTED_ASSOCIATION_ELEMENT_SIZES.contains(keysAssociation.size())) {
                LOGGER.warn("Found " + keysAssociation.size() + " elements in a vicc entry rather than the expected "
                        + EXPECTED_ASSOCIATION_ELEMENT_SIZES);
                LOGGER.warn(keysAssociation);
            } else {
                viccEntryBuilder.association(createAssociation(viccEntryElement));
            }

            viccEntryBuilder.tags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("tags")));
            viccEntryBuilder.devTags(jsonArrayToStringList(viccEntryElement.getAsJsonArray("dev_tags")));

            entries.add(viccEntryBuilder.build());
        }
        reader.close();

        return entries;
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
                .symbol(objectGeneIdentifier.getAsJsonPrimitive("symbol").toString())
                .entrezId(objectGeneIdentifier.getAsJsonPrimitive("entrez_id").toString())
                .ensemblGeneId(objectGeneIdentifier.getAsJsonPrimitive("ensembl_gene_id").toString())
                .build();
    }

    @NotNull
    private static Association createAssociation(JsonObject viccEntryElement) {
        JsonElement elementAssociation = viccEntryElement.get("association");
        LOGGER.info(elementAssociation.getAsJsonObject().keySet());
        Association associationBuilder = ImmutableAssociation.builder()
                .variantName(Strings.EMPTY)
                .evidence(Lists.newArrayList())
                .evidenceLevel(Strings.EMPTY)
                .evidenceLabel(Strings.EMPTY)
                .responseType(Strings.EMPTY)
                .drugLabels(Strings.EMPTY)
                .sourceLink(Strings.EMPTY)
                .publicationUrls(Lists.newArrayList())
                .phenotype(ImmutablePhenotype.builder()
                        .type(ImmutablePhenotypeType.builder().source(Strings.EMPTY).term(Strings.EMPTY).id(Strings.EMPTY).build())
                        .description(Strings.EMPTY)
                        .family(Strings.EMPTY)
                        .build())
                .description(Strings.EMPTY)
                .environmentalContexts(Lists.newArrayList())
                .oncogenic(Strings.EMPTY)
                .build();

        for (String keysAssociation : elementAssociation.getAsJsonObject().keySet()) {
            LOGGER.info(elementAssociation.getAsJsonObject().keySet());
            associationBuilder = ImmutableAssociation.builder()
                    .variantName(elementAssociation.getAsJsonObject().has("variant_name") ? elementAssociation.getAsJsonObject()
                            .getAsJsonPrimitive("variant_name")
                            .toString() : null)

                    .evidence(Lists.newArrayList(createEvidence(elementAssociation)))
                    .evidenceLevel(elementAssociation.getAsJsonObject().getAsJsonPrimitive("evidence_level").toString())
                    .evidenceLabel(elementAssociation.getAsJsonObject().getAsJsonPrimitive("evidence_label").toString())
                    .responseType(elementAssociation.getAsJsonObject().getAsJsonPrimitive("response_type").toString())
                    .drugLabels(elementAssociation.getAsJsonObject().getAsJsonPrimitive("drug_labels").toString())
                    .sourceLink(elementAssociation.getAsJsonObject().has("source_link") ? elementAssociation.getAsJsonObject()
                            .getAsJsonPrimitive("source_link")
                            .toString() : null)
                    .publicationUrls(Lists.newArrayList(elementAssociation.getAsJsonObject()
                            .getAsJsonPrimitive("publication_url")
                            .toString()))
                    .phenotype(createPhenotype(elementAssociation))
                    .description(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString())
                    .environmentalContexts(Lists.newArrayList(createEnvironmentalContexts(elementAssociation)))
                    .oncogenic(elementAssociation.getAsJsonObject().has("oncogenic") ? elementAssociation.getAsJsonObject()
                            .getAsJsonPrimitive("oncogenic")
                            .toString() : null)
                    .build();
        }

        return associationBuilder;
    }

    private static EnvironmentalContext createEnvironmentalContexts(JsonElement elementAssociation) {
        ImmutableEnvironmentalContext builderEnvironmentalContext = ImmutableEnvironmentalContext.builder()
                .term(Strings.EMPTY)
                .description(Strings.EMPTY)
                .taxonomy(ImmutableTaxonomy.builder()
                        .kingdom(Strings.EMPTY)
                        .directParent(Strings.EMPTY)
                        .classs(Strings.EMPTY)
                        .subClass(Strings.EMPTY)
                        .superClass(Strings.EMPTY)
                        .build())
                .source(Strings.EMPTY)
                .usanStem(Strings.EMPTY)
                .approvedCountries(Lists.newArrayList())
                .id(Strings.EMPTY)
                .build();

        for (JsonElement elementEnvironmentContext : elementAssociation.getAsJsonObject().get("environmentalContexts").getAsJsonArray()) {
            StringBuilder approvedCountriesElementString = new StringBuilder();
            String approvedCountriesAdding = Strings.EMPTY;
            if (elementEnvironmentContext.getAsJsonObject().has("approved_countries")) {
                for (JsonElement approvedCountriesElement : elementEnvironmentContext.getAsJsonObject()
                        .get("approved_countries")
                        .getAsJsonArray()) {
                    approvedCountriesElementString.append(approvedCountriesElement).append(",");
                    approvedCountriesAdding = approvedCountriesElementString.substring(0, approvedCountriesElementString.length()-1);
                }
            }

            builderEnvironmentalContext = ImmutableEnvironmentalContext.builder()
                    .term(elementEnvironmentContext.getAsJsonObject().has("term") ? elementEnvironmentContext.getAsJsonObject()
                            .get("term")
                            .getAsString() : null)
                    .description(elementEnvironmentContext.getAsJsonObject().get("description").getAsString())
                    .taxonomy(elementEnvironmentContext.getAsJsonObject().has("taxonomy")
                            ? createTaxonomy(elementEnvironmentContext.getAsJsonObject().get("taxonomy"))
                            : null)
                    .source(elementEnvironmentContext.getAsJsonObject().has("source") ? elementEnvironmentContext.getAsJsonObject()
                            .get("source")
                            .getAsString() : null)
                    .usanStem(elementEnvironmentContext.getAsJsonObject().has("usan_stem") ? elementEnvironmentContext.getAsJsonObject()
                            .get("usan_stem")
                            .getAsString() : null)
                    .approvedCountries(elementEnvironmentContext.getAsJsonObject().has("approved_countries") ? Lists.newArrayList(
                            approvedCountriesAdding) : null)
                    .id(elementEnvironmentContext.getAsJsonObject().has("id") ? elementEnvironmentContext.getAsJsonObject()
                            .get("id")
                            .getAsString() : null)
                    .build();
        } return builderEnvironmentalContext;
    }

    private static Taxonomy createTaxonomy(JsonElement elementEnvironmentContext) {
        return ImmutableTaxonomy.builder()
                .kingdom(elementEnvironmentContext.getAsJsonObject().get("kingdom").getAsString())
                .directParent(elementEnvironmentContext.getAsJsonObject().get("direct-parent").getAsString())
                .classs(elementEnvironmentContext.getAsJsonObject().get("class").getAsString())
                .subClass(elementEnvironmentContext.getAsJsonObject().has("subclass") ? elementEnvironmentContext.getAsJsonObject()
                        .get("subclass")
                        .getAsString() : null)
                .superClass(elementEnvironmentContext.getAsJsonObject().get("superclass").getAsString())
                .build();
    }

    private static Evidence createEvidence(JsonElement elementAssociation) {
        ImmutableEvidence builderEvidence = ImmutableEvidence.builder()
                .info(ImmutableEvidenceInfo.builder().publications(Lists.newArrayList()).build())
                .evidenceType(ImmutableEvidenceType.builder().id(Strings.EMPTY).sourceName(Strings.EMPTY).build())
                .description(Strings.EMPTY)
                .build();

        for (JsonElement elementEvidence : elementAssociation.getAsJsonObject().get("evidence").getAsJsonArray()) {
            builderEvidence = ImmutableEvidence.builder()
                    .info(createEvidenceInfo(elementEvidence.getAsJsonObject().get("info")))
                    .evidenceType(createEvidenceType(elementEvidence.getAsJsonObject().get("evidenceType")))
                    .description(elementEvidence.getAsJsonObject().get("description").getAsString())
                    .build();
        }
        return builderEvidence;
    }

    private static EvidenceType createEvidenceType(JsonElement elementAssociation) {
        return ImmutableEvidenceType.builder()
                .sourceName(elementAssociation.getAsJsonObject().get("sourceName").toString())
                .id(elementAssociation.getAsJsonObject().has("id") ? elementAssociation.getAsJsonObject().get("id").toString() : null)
                .build();
    }

    private static EvidenceInfo createEvidenceInfo(JsonElement elementAssociation) {
        return ImmutableEvidenceInfo.builder()
                .publications(Lists.newArrayList(elementAssociation.getAsJsonObject().get("publications").toString()))
                .build();
    }

    private static PhenotypeType createPhenotypeType(JsonElement elementAssociation) {
        return ImmutablePhenotypeType.builder()
                .source(elementAssociation.getAsJsonObject().get("source").getAsString())
                .term(elementAssociation.getAsJsonObject().get("term").getAsString())
                .id(elementAssociation.getAsJsonObject().get("id").getAsString())
                .build();
    }

    private static Phenotype createPhenotype(JsonElement elementAssociation) {
        return ImmutablePhenotype.builder()
                .type(createPhenotypeType(elementAssociation.getAsJsonObject().get("phenotype").getAsJsonObject().get("type")))
                .description(elementAssociation.getAsJsonObject().get("phenotype").getAsJsonObject().get("description").toString())
                .family(elementAssociation.getAsJsonObject().get("phenotype").getAsJsonObject().get("family").toString())
                .build();
    }
}
