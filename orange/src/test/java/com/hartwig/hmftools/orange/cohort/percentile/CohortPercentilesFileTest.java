package com.hartwig.hmftools.orange.cohort.percentile;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Collection;

import com.google.common.collect.Multimap;
import com.google.common.io.Resources;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CohortPercentilesFileTest {

    private static final String EXAMPLE_PERCENTILE_TSV = Resources.getResource("cohort/percentile/example_percentile_file.tsv").getPath();

    @Test
    public void canReadExampleCohortPercentilesTsv() throws IOException {
        Multimap<PercentileType, CohortPercentiles> map = CohortPercentilesFile.read(EXAMPLE_PERCENTILE_TSV);

        assertEquals(1, map.keySet().size());
        assertEquals(2, map.values().size());

        CohortPercentiles ovary = find(map.get(PercentileType.SV_PASS_COUNT), "Ovary");
        assertEquals(5, ovary.percentileValues().size());
    }

    @NotNull
    private static CohortPercentiles find(@NotNull Collection<CohortPercentiles> percentiles, @NotNull String cancerType) {
        for (CohortPercentiles percentile : percentiles) {
            if (percentile.cancerType().equals(cancerType)) {
                return percentile;
            }
        }

        throw new IllegalStateException("Could not find percentile for cancer type: " + cancerType);
    }
}