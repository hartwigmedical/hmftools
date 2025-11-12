package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.interpretation.Drivers;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.orange.algo.util.PurpleDriverTestFactory;

import org.junit.Test;

public class DriversTest
{
    @Test
    public void canSelectNonCanonicalMutationEntries()
    {
        PurpleDriver canonicalMutation = PurpleDriverTestFactory.builder().type(PurpleDriverType.MUTATION).isCanonical(true).build();
        PurpleDriver nonCanonicalMutation = PurpleDriverTestFactory.builder().type(PurpleDriverType.MUTATION).isCanonical(false).build();
        PurpleDriver nonCanonicalAmp = PurpleDriverTestFactory.builder().type(PurpleDriverType.AMP).isCanonical(false).build();

        List<PurpleDriver> drivers = Lists.newArrayList(canonicalMutation, nonCanonicalMutation, nonCanonicalAmp);
        List<PurpleDriver> nonCanonicalMutationDrivers = Drivers.nonCanonicalMutationEntries(drivers);

        assertEquals(1, nonCanonicalMutationDrivers.size());
        assertTrue(nonCanonicalMutationDrivers.contains(nonCanonicalMutation));
    }

    @Test
    public void canSelectCanonicalMutationEntryForGene()
    {
        PurpleDriver canonicalMatchLowDL = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.3)
                .build();

        PurpleDriver canonicalMatchHighDL = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.4)
                .build();

        PurpleDriver nonCanonicalMatch = PurpleDriverTestFactory.builder()
                .type(PurpleDriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(false)
                .driverLikelihood(0.5)
                .build();

        PurpleDriver canonicalOtherDriver =
                PurpleDriverTestFactory.builder()
                        .type(PurpleDriverType.AMP)
                        .gene("gene 1")
                        .isCanonical(true)
                        .driverLikelihood(0.6)
                        .build();

        PurpleDriver canonicalOtherGene =
                PurpleDriverTestFactory.builder()
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