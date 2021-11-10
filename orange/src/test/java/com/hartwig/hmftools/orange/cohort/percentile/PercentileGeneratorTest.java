package com.hartwig.hmftools.orange.cohort.percentile;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableObservation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class PercentileGeneratorTest {

    private static final double EPSILON = 1.0E-10;

    @Test
    public void canGeneratePercentilesForEmptyData() {
        PercentileGenerator generator = new PercentileGenerator(sample -> "type");

        assertTrue(generator.run(Lists.newArrayList()).isEmpty());
    }

    @Test
    public void canGeneratePercentilesForLimitedData() {
        PercentileGenerator generator = new PercentileGenerator(sample -> "type");

        ImmutableObservation.Builder builder = createTestBuilder();

        List<Observation> observations = Lists.newArrayList();
        observations.add(builder.value(1).build());
        observations.add(builder.value(2).build());
        observations.add(builder.value(3).build());
        observations.add(builder.value(4).build());
        observations.add(builder.value(5).build());

        List<CohortPercentiles> percentiles = generator.run(observations);
        assertEquals(2, percentiles.size());
        assertEquals(5, percentiles.get(0).cohortSize());

        assertPercentilesForObservations(observations, percentiles.get(0).values());
    }

    @Test
    public void doNotCrashOnMissingCancerTypes() {
        PercentileGenerator generator = new PercentileGenerator(sample -> null);

        ImmutableObservation.Builder builder = createTestBuilder();

        List<Observation> observations = Lists.newArrayList(builder.value(1).build());

        assertNotNull(generator.run(observations));
    }

    @Test
    public void canGeneratePercentilesForExactBucketCount() {
        PercentileGenerator generator = new PercentileGenerator(sample -> "type");

        ImmutableObservation.Builder builder = createTestBuilder();

        List<Observation> observations = Lists.newArrayList();
        for (int i = 0; i < PercentileConstants.BUCKET_COUNT; i++) {
            observations.add(builder.value(i+1).build());
        }

        List<CohortPercentiles> percentiles = generator.run(observations);
        assertEquals(2, percentiles.size());
        assertEquals(PercentileConstants.BUCKET_COUNT, percentiles.get(0).cohortSize());

        assertPercentilesForObservations(observations, percentiles.get(0).values());
    }

    @Test
    public void canGeneratePercentilesForMoreDataThanBucketCount() {
        PercentileGenerator generator = new PercentileGenerator(sample -> "type");

        ImmutableObservation.Builder builder = createTestBuilder();

        List<Observation> observations = Lists.newArrayList();
        int count = PercentileConstants.BUCKET_COUNT * 20;
        for (int i = 0; i < count; i++) {
            observations.add(builder.value(i+1).build());
        }

        List<CohortPercentiles> percentiles = generator.run(observations);
        assertEquals(2, percentiles.size());
        assertEquals(count, percentiles.get(0).cohortSize());

        assertPercentilesForObservations(observations, percentiles.get(0).values());
    }

    @NotNull
    private static ImmutableObservation.Builder createTestBuilder() {
        return ImmutableObservation.builder().type(PercentileType.SV_TMB).sample(ImmutableSample.builder().sampleId("test").build());
    }

    private static void assertPercentilesForObservations(@NotNull List<Observation> observations, @NotNull List<Double> percentiles) {
        // Percentiles should always have the same number of buckets.
        assertEquals(PercentileConstants.BUCKET_COUNT, percentiles.size());

        // The first element of percentiles should always match the lowest number.
        assertEquals(observations.get(0).value(), percentiles.get(0), EPSILON);

        // The last element of percentiles should always match the highest number.
        assertEquals(observations.get(observations.size() - 1).value(), percentiles.get(percentiles.size() - 1), EPSILON);
    }
}