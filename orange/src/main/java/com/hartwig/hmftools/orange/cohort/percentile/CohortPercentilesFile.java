package com.hartwig.hmftools.orange.cohort.percentile;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.createFieldsIndexMap;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;

import org.jetbrains.annotations.NotNull;

public final class CohortPercentilesFile {

    private static final String LINE_DELIMITER = "\t";
    private static final String PERCENTILE_DELIMITER = ";";

    private CohortPercentilesFile() {
    }

    @NotNull
    public static Multimap<PercentileType, CohortPercentiles> read(@NotNull String tsv) throws IOException {
        Multimap<PercentileType, CohortPercentiles> map = ArrayListMultimap.create();

        List<String> lines = Files.readAllLines(new File(tsv).toPath());

        String header = lines.get(0);
        Map<String, Integer> fields = createFieldsIndexMap(header, LINE_DELIMITER);
        lines.remove(0);

        for (String line : lines) {
            String[] values = line.split(LINE_DELIMITER, -1);

            CohortPercentiles percentiles = ImmutableCohortPercentiles.builder()
                    .cancerType(values[fields.get("cancerType")])
                    .cohortSize(Integer.parseInt(values[fields.get("cohortSize")]))
                    .percentileValues(toPercentileValues(values[fields.get("percentiles")]))
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
