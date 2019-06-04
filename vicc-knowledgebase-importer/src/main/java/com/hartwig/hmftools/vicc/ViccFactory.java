package com.hartwig.hmftools.vicc;

import java.io.FileReader;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
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

    private ViccFactory() {
    }

    public static List<ViccEntry> readViccKnowledgebaseJsonFile(@NotNull String jsonPath) throws IOException {
        JsonParser parser = new JsonParser();
        JsonReader reader = new JsonReader(new FileReader(jsonPath));
        reader.setLenient(true);

        List<ViccEntry> entries = Lists.newArrayList();
        while (reader.peek() != JsonToken.END_DOCUMENT) {
            JsonObject object = parser.parse(reader).getAsJsonObject();

            LOGGER.info("Found " + object.size() + " elements");

            ImmutableViccEntry.Builder viccEntryBuilder = createPreGeneratedBuilder();
            viccEntryBuilder.source(object.getAsJsonPrimitive("source").getAsString());
            entries.add(viccEntryBuilder.build());
        }
        reader.close();

        return entries;
    }

    @NotNull
    private static ImmutableViccEntry.Builder createPreGeneratedBuilder() {
        ImmutableViccEntry.Builder viccEntryBuilder = ImmutableViccEntry.builder();
        viccEntryBuilder.source(Strings.EMPTY);
        viccEntryBuilder.genes(Lists.newArrayList());
        viccEntryBuilder.geneIdentifiers(Lists.newArrayList());
        viccEntryBuilder.featureNames(Lists.newArrayList());
        viccEntryBuilder.features(Lists.newArrayList());
        viccEntryBuilder.association(createAssociation());
        viccEntryBuilder.tags(Lists.newArrayList());
        viccEntryBuilder.devTags(Lists.newArrayList());
        return viccEntryBuilder;
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
