package com.hartwig.hmftools.vicc;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

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
            viccEntryBuilder.association(createAssociation(viccEntryElement));

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
        JsonElement geneIdentifiers = viccEntryElement.get("gene_identifiers");
        List<GeneIdentifier> listGeneIdentifiers = Lists.newArrayList();

        for (JsonElement elementGeneIdentifier : geneIdentifiers.getAsJsonArray()) {
            listGeneIdentifiers.add(toGeneIdentifier(elementGeneIdentifier));

        }
        return listGeneIdentifiers;
    }

    @NotNull
    private static GeneIdentifier toGeneIdentifier(@NotNull JsonElement elementGeneIdentifier) {
        return ImmutableGeneIdentifier.builder()
                .symbol(elementGeneIdentifier.getAsJsonObject().getAsJsonPrimitive("symbol").toString())
                .entrezId(elementGeneIdentifier.getAsJsonObject().getAsJsonPrimitive("entrez_id").toString())
                .ensemblGeneId(elementGeneIdentifier.getAsJsonObject().getAsJsonPrimitive("ensembl_gene_id").toString())
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
                        .id(Strings.EMPTY)
                        .build())
                .description(Strings.EMPTY)
                .environmentalContexts(Lists.newArrayList())
                .build();

        for (String keysAssociation : elementAssociation.getAsJsonObject().keySet()) {
            LOGGER.info(elementAssociation.getAsJsonObject().keySet());
            associationBuilder = ImmutableAssociation.builder()
                    .variantName(Strings.EMPTY)
                    .evidence(Lists.newArrayList())
                    .evidenceLevel(elementAssociation.getAsJsonObject().getAsJsonPrimitive("evidence_level").toString())
                    .evidenceLabel(elementAssociation.getAsJsonObject().getAsJsonPrimitive("evidence_label").toString())
                    .responseType(elementAssociation.getAsJsonObject().getAsJsonPrimitive("response_type").toString())
                    .drugLabels(elementAssociation.getAsJsonObject().getAsJsonPrimitive("drug_labels").toString())
                    .sourceLink(Strings.EMPTY)
                    .publicationUrls(Lists.newArrayList(elementAssociation.getAsJsonObject()
                            .getAsJsonPrimitive("publication_url")
                            .toString()))
                    .phenotype(createPhenotype(elementAssociation))
                    .description(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString())
                    .environmentalContexts(Lists.newArrayList(createEnvironmentalContexts(elementAssociation)))
                    .build();
        }

        return associationBuilder;
    }

    private static EnvironmentalContext createEnvironmentalContexts(JsonElement elementAssociation) {
        return ImmutableEnvironmentalContext.builder()
                .term(Strings.EMPTY)
                .description(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString())
                .taxonomy(createTaxonomy(elementAssociation))
                .source(Strings.EMPTY)
                .usanStem(Strings.EMPTY)
                .approvedCountries(Lists.newArrayList())
                .id(Strings.EMPTY)
                .build();
    }

    private static Taxonomy createTaxonomy(JsonElement elementAssociation) {
        return ImmutableTaxonomy.builder()
                .kingdom(Strings.EMPTY)
                .directParent(Strings.EMPTY)
                .classs(Strings.EMPTY)
                .subClass(Strings.EMPTY)
                .superClass(Strings.EMPTY)
                .build();
    }

    private static Evidence createEvidence(JsonElement elementAssociation) {
        LOGGER.info(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString());
        return ImmutableEvidence.builder()
                .info(createEvidenceInfo(elementAssociation))
                .evidenceType(createEvidenceType(elementAssociation))
                .description(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString())
                .build();
    }

    private static EvidenceType createEvidenceType(JsonElement elementAssociation) {
        return ImmutableEvidenceType.builder()
                .sourceName(elementAssociation.getAsJsonObject().getAsJsonPrimitive("sourceName").toString())
                .id(elementAssociation.getAsJsonObject().getAsJsonPrimitive("id").toString())
                .build();
    }

    private static EvidenceInfo createEvidenceInfo(JsonElement elementAssociation) {
        return ImmutableEvidenceInfo.builder()
                .publications(Lists.newArrayList(elementAssociation.getAsJsonObject().getAsJsonPrimitive("publications").toString()))
                .build();
    }

    private static Phenotype createPhenotype(JsonElement elementAssociation) {
        return ImmutablePhenotype.builder()
                .type(ImmutablePhenotypeType.builder()
                        .source(Strings.EMPTY)
                        .term(Strings.EMPTY)
                        .id(Strings.EMPTY)
                        .build())
                .description(elementAssociation.getAsJsonObject().getAsJsonPrimitive("description").toString())
                .family(Strings.EMPTY)
                .id(Strings.EMPTY)
                .build();
    }
}
