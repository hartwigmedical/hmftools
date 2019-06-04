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
import com.hartwig.hmftools.vicc.datamodel.ImmutableAssociation;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotype;
import com.hartwig.hmftools.vicc.datamodel.ImmutablePhenotypeType;
import com.hartwig.hmftools.vicc.datamodel.ImmutableViccEntry;
import com.hartwig.hmftools.vicc.datamodel.Phenotype;
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

            viccEntryBuilder.geneIdentifiers(Lists.newArrayList());

            if (viccEntryElement.has("feature_names")) {
                JsonElement featureNames = viccEntryElement.get("feature_names");
                if (featureNames.isJsonArray()) {
                    viccEntryBuilder.featureNames(jsonArrayToStringList(featureNames.getAsJsonArray()));
                } else if (featureNames.isJsonPrimitive()) {
                    viccEntryBuilder.featureNames(Lists.newArrayList(featureNames.getAsJsonPrimitive().getAsString()));
                }
            }

            viccEntryBuilder.features(Lists.newArrayList());
            viccEntryBuilder.association(createAssociation());

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
    private static Association createAssociation() {
        return ImmutableAssociation.builder()
                .variantName(Strings.EMPTY)
                .evidence(Lists.newArrayList())
                .evidenceLevel(Strings.EMPTY)
                .evidenceLabel(Strings.EMPTY)
                .responseType(Strings.EMPTY)
                .drugLabels(Strings.EMPTY)
                .sourceLink(Strings.EMPTY)
                .publicationUrls(Lists.newArrayList())
                .phenotype(createPhenotype())
                .description(Strings.EMPTY)
                .environmentalContexts(Lists.newArrayList())
                .build();
    }

    private static Phenotype createPhenotype() {
        return ImmutablePhenotype.builder()
                .type(ImmutablePhenotypeType.builder().source(Strings.EMPTY).term(Strings.EMPTY).id(Strings.EMPTY).build())
                .description(Strings.EMPTY)
                .family(Strings.EMPTY)
                .id(Strings.EMPTY)
                .build();
    }
}
