package com.hartwig.hmftools.orange.cohort.percentile;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.StringJoiner;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public final class CohortPercentilesFile {

    private static final String COHORT_FILE_NAME = "orange_cohort_percentiles.tsv";

    private static final String FIELD_DELIMITER = "\t";
    private static final String PERCENTILE_DELIMITER = ";";

    private CohortPercentilesFile() {
    }

    @NotNull
    public static String generateOutputTsv(@NotNull String outputDirectory) {
        String path = outputDirectory.endsWith(File.separator) ? outputDirectory : outputDirectory + File.separator;
        return path + COHORT_FILE_NAME;
    }

    public static void write(@NotNull String tsv, @NotNull Multimap<PercentileType, CohortPercentiles> percentileMap) throws IOException {
        List<String> lines = Lists.newArrayList();
        lines.add(header());
        lines.addAll(toLines(percentileMap));

        Files.write(new File(tsv).toPath(), lines);
    }

    @NotNull
    public static Multimap<PercentileType, CohortPercentiles> read(@NotNull String tsv) throws IOException {
        List<String> lines = Files.readAllLines(new File(tsv).toPath());

        Map<String, Integer> fields = createFieldsIndexMap(lines.get(0), FIELD_DELIMITER);

        return fromLines(fields, lines.subList(1, lines.size()));
    }

    @NotNull
    @VisibleForTesting
    static String header() {
        return new StringJoiner(FIELD_DELIMITER).add("type").add("cancerType").add("cohortSize").add("percentiles").toString();
    }

    @NotNull
    @VisibleForTesting
    static List<String> toLines(@NotNull Multimap<PercentileType, CohortPercentiles> percentileMap) {
        List<String> lines = Lists.newArrayList();
        for (Map.Entry<PercentileType, Collection<CohortPercentiles>> entry : percentileMap.asMap().entrySet()) {
            for (CohortPercentiles percentiles : entry.getValue()) {
                lines.add(toLine(entry.getKey(), percentiles));
            }
        }
        return lines;
    }

    @NotNull
    private static String toLine(@NotNull PercentileType type, @NotNull CohortPercentiles percentiles) {
        StringJoiner percentileField = new StringJoiner(PERCENTILE_DELIMITER);
        for (double value : percentiles.values()) {
            percentileField.add(String.valueOf(value));
        }

        return new StringJoiner(FIELD_DELIMITER).add(type.toString())
                .add(percentiles.cancerType())
                .add(String.valueOf(percentiles.cohortSize()))
                .add(percentileField.toString())
                .toString();
    }

    @NotNull
    static Multimap<PercentileType, CohortPercentiles> fromLines(@NotNull Map<String, Integer> fields, @NotNull List<String> lines) {
        Multimap<PercentileType, CohortPercentiles> map = ArrayListMultimap.create();

        for (String line : lines) {
            String[] values = line.split(FIELD_DELIMITER, -1);

            CohortPercentiles percentiles = ImmutableCohortPercentiles.builder()
                    .cancerType(values[fields.get("cancerType")])
                    .cohortSize(Integer.parseInt(values[fields.get("cohortSize")]))
                    .values(toPercentileValues(values[fields.get("percentiles")]))
                    .build();

            PercentileType type = PercentileType.valueOf(values[fields.get("type")]);
            if (containsCancerTypeForPercentileType(map.get(type), percentiles.cancerType())) {
                throw new IllegalStateException("Double entry of cancer type " + percentiles.cancerType() + " for percentile type " + type);
            }

            map.put(type, percentiles);
        }
        return map;
    }

    private static boolean containsCancerTypeForPercentileType(@NotNull Collection<CohortPercentiles> percentiles,
            @NotNull String cancerType) {
        for (CohortPercentiles percentile : percentiles) {
            if (percentile.cancerType().equals(cancerType)) {
                return true;
            }
        }

        return false;
    }

    @NotNull
    private static List<Double> toPercentileValues(@NotNull String percentileString) {
        List<Double> percentileValues = Lists.newArrayList();
        for (String percentile : percentileString.split(PERCENTILE_DELIMITER)) {
            percentileValues.add(Double.parseDouble(percentile));
        }
        return percentileValues;
    }
}
