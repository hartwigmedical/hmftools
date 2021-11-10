package com.hartwig.hmftools.orange.cohort.percentile;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableObservation;
import com.hartwig.hmftools.orange.cohort.datamodel.ImmutableSample;
import com.hartwig.hmftools.orange.cohort.datamodel.Observation;

import org.junit.Test;

public class PercentileGeneratorTest {

    @Test
    public void canGeneratePercentiles() {
        PercentileGenerator generator = new PercentileGenerator(sample -> "type");

        assertTrue(generator.run(Lists.newArrayList()).isEmpty());

        ImmutableObservation.Builder builder =
                ImmutableObservation.builder().type(PercentileType.SV_TMB).sample(ImmutableSample.builder().sampleId("test").build());

        List<Observation> observationLimited = Lists.newArrayList();
        observationLimited.add(builder.value(1).build());
        observationLimited.add(builder.value(2).build());
        observationLimited.add(builder.value(3).build());
        observationLimited.add(builder.value(4).build());
        observationLimited.add(builder.value(5).build());

        List<CohortPercentiles> percentilesLimited = generator.run(observationLimited);
        assertEquals(2, percentilesLimited.size());
        assertEquals(5, percentilesLimited.get(0).cohortSize());
        assertEquals(PercentileConstants.BUCKET_COUNT, percentilesLimited.get(0).values().size());

        List<Observation> observationsMassive = Lists.newArrayList();
        int massiveCount = PercentileConstants.BUCKET_COUNT * 3;
        for (int i = 0; i < massiveCount; i++) {
            observationsMassive.add(builder.value(i).build());
        }

        List<CohortPercentiles> percentilesMassive = generator.run(observationsMassive);
        assertEquals(2, percentilesMassive.size());
        assertEquals(massiveCount, percentilesMassive.get(0).cohortSize());
        assertEquals(PercentileConstants.BUCKET_COUNT, percentilesMassive.get(0).values().size());
    }
}