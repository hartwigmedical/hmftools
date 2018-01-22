package com.hartwig.hmftools.common.cosmic.fusions;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.exception.EmptyFileException;
import com.hartwig.hmftools.common.io.reader.FileReader;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class COSMICGeneFusions {

    private static final Logger LOGGER = LogManager.getLogger(COSMICGeneFusions.class);

    private static final int FIVE_COLUMN = 0;
    private static final int THREE_COLUMN = 1;
    private static final int URL_COLUMN = 2;
    private static final int FIELD_COUNT = 3;

    private static final String FIELD_SEPARATOR = ",";

    private static final int MIN_PROMISCUOUS_PARTNER = 3;

    private COSMICGeneFusions() {
    }

    @NotNull
    public static COSMICGeneFusionModel readFromCSV(@NotNull final String pathToCSV) throws IOException, EmptyFileException {
        final List<COSMICGeneFusionData> items = Lists.newArrayList();
        final List<String> lines = FileReader.build().readLines(new File(pathToCSV).toPath());
        lines.remove(0); // NERA: delete header

        final Map<String, Integer> countThree = Maps.newHashMap();
        final Map<String, Integer> countFive = Maps.newHashMap();

        for (final String line : lines) {
            final String[] parts = line.split(FIELD_SEPARATOR, FIELD_COUNT);
            if (parts.length == FIELD_COUNT) {
                items.add(ImmutableCOSMICGeneFusionData.of(gene(parts[FIVE_COLUMN]),
                        transcript(parts[FIVE_COLUMN]),
                        gene(parts[THREE_COLUMN]),
                        transcript(parts[THREE_COLUMN]),
                        parts[URL_COLUMN]));
                countFive.compute(gene(parts[FIVE_COLUMN]), (k, v) -> v == null ? 1 : v + 1);
                countThree.compute(gene(parts[THREE_COLUMN]), (k, v) -> v == null ? 1 : v + 1);
            } else if (parts.length > 0) {
                LOGGER.warn("Could not properly parse line in gene fusion csv: " + line);
            }
        }

        return ImmutableCOSMICGeneFusionModel.of(items,
                countFive.entrySet()
                        .stream()
                        .filter(e -> e.getValue() > MIN_PROMISCUOUS_PARTNER)
                        .map(e -> ImmutableCOSMICPromiscuousGene.of(gene(e.getKey()), transcript(e.getKey())))
                        .collect(Collectors.toList()),
                countThree.entrySet()
                        .stream()
                        .filter(e -> e.getValue() > MIN_PROMISCUOUS_PARTNER)
                        .map(e -> ImmutableCOSMICPromiscuousGene.of(gene(e.getKey()), transcript(e.getKey())))
                        .collect(Collectors.toList()));
    }

    @NotNull
    private static String gene(@NotNull String input) {
        return input.contains("_") ? input.split("_")[0] : input;
    }

    @Nullable
    private static String transcript(@NotNull String input) {
        if (input.contains("_")) {
            final String transcript = input.split("_", 2)[1];
            // NERA: only want ENSEMBL transcript ID
            return transcript.startsWith("ENST") ? transcript : null;
        }
        return null;
    }
}
