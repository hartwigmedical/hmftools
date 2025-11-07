package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.drivers.DndsCalculator.probabilityDriverVariant;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.Collections;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCatalogFactory;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.purple.DriverGeneResource;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Ignore;
import org.junit.Test;

public class DndsCalculatorTest
{
    private static final double EPSILON = 0.0001;

    @Ignore
    @Test
    public void testHIST2H3DMissense()
    {
        DndsDriverImpactLikelihood driverImpactLikelihood = buildImpactLikelihood();
        double value = probabilityDriverVariant(27742, driverImpactLikelihood); // mDriverGenes.OncoLikelihoodMap.get("H3C13"));
        assertEquals(0.6109, value, EPSILON);
    }

    @Ignore
    @Test
    public void testABL1Missense()
    {
        DndsDriverImpactLikelihood driverImpactLikelihood = buildImpactLikelihood(); // onco.get("ABL1")
        double value = probabilityDriverVariant(996698, driverImpactLikelihood);
        assertEquals(0.0039, value, EPSILON);
    }

    @Ignore
    @Test
    public void testGATA3Indel()
    {
        DndsDriverImpactLikelihood driverImpactLikelihood = buildImpactLikelihood(); // tsg.get("GATA3").indel())
        double value = probabilityDriverVariant(587, driverImpactLikelihood);
        assertEquals(0.9956, value, EPSILON);
    }

    @Test
    public void testMultipleZeroNonDriver()
    {
        DndsDriverImpactLikelihood indelLikelihood = testBuilder().driversPerSample(0.01).passengersPerMutation(0).build();
        double value = probabilityDriverVariant(1000, 1000, indelLikelihood, indelLikelihood);
        assertEquals(0, value, EPSILON);
    }

    @Test
    public void testFallBackOnSingleProbabilityIfMultiFailsDueToZeroValues()
    {
        DndsDriverImpactLikelihood nonsense = testBuilder().driversPerSample(2E-4).passengersPerMutation(4E-9).build();
        double singleNonsenseLikelihood = probabilityDriverVariant(1, nonsense);
        assertTrue(singleNonsenseLikelihood > 0.1);

        DndsDriverImpactLikelihood splice = testBuilder().driversPerSample(0).passengersPerMutation(0).build();
        assertEquals(0, probabilityDriverVariant(1, splice), EPSILON);

        double victim = probabilityDriverVariant(1, 1, nonsense, splice);
        assertEquals(singleNonsenseLikelihood, victim, EPSILON);

        victim = probabilityDriverVariant(1, 1, splice, nonsense);
        assertEquals(singleNonsenseLikelihood, victim, EPSILON);
    }

    @Test
    public void testZeroNonDriverWithStandard()
    {
        DndsDriverImpactLikelihood indelLikelihood = testBuilder().driversPerSample(0.01).passengersPerMutation(0).build();
        DndsDriverImpactLikelihood missenseLikelihood = testBuilder().driversPerSample(0.01).passengersPerMutation(10E-8).build();

        double expectedMissense = probabilityDriverVariant(10000, missenseLikelihood);
        double value = probabilityDriverVariant(10000, 1000, missenseLikelihood, indelLikelihood);
        assertEquals(expectedMissense, value, EPSILON);
    }

    private static DndsDriverImpactLikelihood buildImpactLikelihood()
    {
        return ImmutableDndsDriverImpactLikelihood.builder().driversPerSample(1).passengersPerMutation(1).build();
    }

    @NotNull
    private static ImmutableDndsDriverImpactLikelihood.Builder testBuilder()
    {
        return ImmutableDndsDriverImpactLikelihood.builder();
    }
}
