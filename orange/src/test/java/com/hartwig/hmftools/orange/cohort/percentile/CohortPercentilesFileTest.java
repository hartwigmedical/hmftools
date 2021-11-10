package com.hartwig.hmftools.orange.cohort.percentile;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.utils.FileWriterUtils;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CohortPercentilesFileTest {

    private static final String EXAMPLE_PERCENTILE_TSV = Resources.getResource("cohort/percentile/example_cohort_percentiles.tsv").getPath();

    @Test
    public void canReadExampleCohortPercentilesTsv() throws IOException {
        Multimap<PercentileType, CohortPercentiles> map = CohortPercentilesFile.read(EXAMPLE_PERCENTILE_TSV);

        assertEquals(1, map.keySet().size());
        assertEquals(4, map.values().size());

        CohortPercentiles ovary = find(map.get(PercentileType.SV_TMB), "Ovary");
        assertEquals(5, ovary.values().size());

        assertNotNull(find(map.get(PercentileType.SV_TMB), CohortConstants.COHORT_PAN_CANCER));
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

    @Test
    public void canConvertBackAndForth() {
        Multimap<PercentileType, CohortPercentiles> testMap = createTestMap();

        List<String> lines = CohortPercentilesFile.toLines(testMap);

        Map<String, Integer> fields = FileWriterUtils.createFieldsIndexMap(CohortPercentilesFile.header(), "\t");

        Multimap<PercentileType, CohortPercentiles> recreatedMap = CohortPercentilesFile.fromLines(fields, lines);

        assertEquals(testMap, recreatedMap);
    }

    @NotNull
    private static Multimap<PercentileType, CohortPercentiles> createTestMap() {
        Multimap<PercentileType, CohortPercentiles> map = ArrayListMultimap.create();
        map.put(PercentileType.SV_TMB,
                ImmutableCohortPercentiles.builder().cancerType("type 1").cohortSize(12).addValues(1, 2, 3, 4).build());
        map.put(PercentileType.SV_TMB,
                ImmutableCohortPercentiles.builder().cancerType("type 2").cohortSize(12).addValues(1, 2, 3, 4).build());
        return map;
    }
}