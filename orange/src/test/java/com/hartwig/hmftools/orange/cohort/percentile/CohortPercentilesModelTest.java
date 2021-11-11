package com.hartwig.hmftools.orange.cohort.percentile;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableObservation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;
import com.hartwig.hmftools.orange.cohort.mapping.CohortConstants;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CohortPercentilesModelTest {

    private static final String CANCER_TYPE = "type";
    private static final PercentileType PERCENTILE_TYPE = PercentileType.SV_TMB;

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canDeterminePercentileForCohort() {
        List<Double> percentiles = Lists.newArrayList(1D, 3D, 5D, 7D, 9D);

        CohortPercentilesModel model = new CohortPercentilesModel(sample -> CANCER_TYPE, createPercentileMap(percentiles));

        // Lower that lowest percentile
        assertEquals(0D, model.percentile(observedValue(0D)).cancerTypePercentile(), EPSILON);

        // In the middle
        assertEquals(0.5, model.percentile(observedValue(5D)).cancerTypePercentile(), EPSILON);

        // Above the middle
        assertEquals(0.6, model.percentile(observedValue(6D)).cancerTypePercentile(), EPSILON);

        // At the top
        assertEquals(0.9, model.percentile(observedValue(9D)).cancerTypePercentile(), EPSILON);

        // Higher than highest percentile
        assertEquals(1D, model.percentile(observedValue(10D)).cancerTypePercentile(), EPSILON);
    }

    @Test (expected = IllegalStateException.class)
    public void crashOnExpectedCohort() {
        List<Double> percentiles = Lists.newArrayList(1D, 3D, 5D, 7D, 9D);
        CohortPercentilesModel model = new CohortPercentilesModel(sample -> "non-existing cancer type", createPercentileMap(percentiles));

        model.percentile(observedValue(0D));

    }

    private static Multimap<PercentileType, CohortPercentiles> createPercentileMap(@NotNull List<Double> percentiles) {
        Multimap<PercentileType, CohortPercentiles> percentileMap = ArrayListMultimap.create();
        percentileMap.put(PERCENTILE_TYPE,
                ImmutableCohortPercentiles.builder().cancerType(CANCER_TYPE).cohortSize(20).values(percentiles).build());

        percentileMap.put(PERCENTILE_TYPE,
                ImmutableCohortPercentiles.builder()
                        .cancerType(CohortConstants.COHORT_PAN_CANCER)
                        .cohortSize(20)
                        .values(percentiles)
                        .build());

        return percentileMap;
    }

    @NotNull
    private static Observation observedValue(double value) {
        return ImmutableObservation.builder()
                .type(PERCENTILE_TYPE)
                .sample(ImmutableSample.builder().sampleId("test").build())
                .value(value)
                .build();
    }
}