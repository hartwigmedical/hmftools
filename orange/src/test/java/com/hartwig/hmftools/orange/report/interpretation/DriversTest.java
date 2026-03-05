package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.orange.algo.purple.PurpleTestFactory;

import org.junit.Test;

public class DriversTest
{
    @Test
    public void canSelectNonCanonicalMutationEntries()
    {
        PurpleDriver canonicalMutation = PurpleTestFactory.purpleDriverBuilder().type(PurpleDriverType.MUTATION).isCanonical(true).build();
        PurpleDriver nonCanonicalMutation = PurpleTestFactory.purpleDriverBuilder().type(PurpleDriverType.MUTATION).isCanonical(false).build();
        PurpleDriver nonCanonicalAmp = PurpleTestFactory.purpleDriverBuilder().type(PurpleDriverType.AMP).isCanonical(false).build();

        List<PurpleDriver> drivers = Lists.newArrayList(canonicalMutation, nonCanonicalMutation, nonCanonicalAmp);
        List<PurpleDriver> nonCanonicalMutationDrivers = Drivers.nonCanonicalMutationEntries(drivers);

        assertEquals(1, nonCanonicalMutationDrivers.size());
        assertTrue(nonCanonicalMutationDrivers.contains(nonCanonicalMutation));
    }

    @Test
    public void canSelectCanonicalMutationEntryForGene()
    {
        PurpleDriver canonicalMatchLowDL = PurpleTestFactory.purpleDriverBuilder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.3)
                .build();

        PurpleDriver canonicalMatchHighDL = PurpleTestFactory.purpleDriverBuilder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.4)
                .build();

        PurpleDriver nonCanonicalMatch = PurpleTestFactory.purpleDriverBuilder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(false)
                .driverLikelihood(0.5)
                .build();

        PurpleDriver canonicalOtherDriver =
                PurpleTestFactory.purpleDriverBuilder()
                        .type(PurpleDriverType.AMP)
                        .gene("gene 1")
                        .isCanonical(true)
                        .driverLikelihood(0.6)
                        .build();

        PurpleDriver canonicalOtherGene =
                PurpleTestFactory.purpleDriverBuilder()
                        .type(PurpleDriverType.AMP)
                        .gene("gene 2")
                        .isCanonical(true)
                        .driverLikelihood(0.7)
                        .build();

        List<PurpleDriver> drivers =
                Lists.newArrayList(canonicalMatchLowDL, canonicalMatchHighDL, nonCanonicalMatch, canonicalOtherDriver, canonicalOtherGene);

        assertEquals(canonicalMatchHighDL, Drivers.canonicalMutationEntryForGene(drivers, "gene 1"));
        assertNull(Drivers.canonicalMutationEntryForGene(drivers, "gene 2"));
    }
}