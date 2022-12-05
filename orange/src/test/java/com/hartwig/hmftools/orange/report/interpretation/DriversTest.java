package com.hartwig.hmftools.orange.report.interpretation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogTestFactory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.junit.Test;

public class DriversTest {

    @Test
    public void canSelectNonCanonicalMutationEntries() {
        DriverCatalog canonicalMutation = DriverCatalogTestFactory.builder().driver(DriverType.MUTATION).isCanonical(true).build();
        DriverCatalog nonCanonicalMutation = DriverCatalogTestFactory.builder().driver(DriverType.MUTATION).isCanonical(false).build();
        DriverCatalog nonCanonicalAmp = DriverCatalogTestFactory.builder().driver(DriverType.AMP).isCanonical(false).build();

        List<DriverCatalog> drivers = Lists.newArrayList(canonicalMutation, nonCanonicalMutation, nonCanonicalAmp);
        List<DriverCatalog> nonCanonicalMutationDrivers = Drivers.nonCanonicalMutationEntries(drivers);

        assertEquals(1, nonCanonicalMutationDrivers.size());
        assertTrue(nonCanonicalMutationDrivers.contains(nonCanonicalMutation));
    }

    @Test
    public void canSelectCanonicalMutationEntryForGene() {
        DriverCatalog canonicalMatchLowDL = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.3)
                .build();

        DriverCatalog canonicalMatchHighDL = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(true)
                .driverLikelihood(0.4)
                .build();

        DriverCatalog nonCanonicalMatch = DriverCatalogTestFactory.builder()
                .driver(DriverType.MUTATION)
                .gene("gene 1")
                .isCanonical(false)
                .driverLikelihood(0.5)
                .build();

        DriverCatalog canonicalOtherDriver =
                DriverCatalogTestFactory.builder().driver(DriverType.AMP).gene("gene 1").isCanonical(true).driverLikelihood(0.6).build();

        DriverCatalog canonicalOtherGene =
                DriverCatalogTestFactory.builder().driver(DriverType.AMP).gene("gene 2").isCanonical(true).driverLikelihood(0.7).build();

        List<DriverCatalog> drivers =
                Lists.newArrayList(canonicalMatchLowDL, canonicalMatchHighDL, nonCanonicalMatch, canonicalOtherDriver, canonicalOtherGene);

        assertEquals(canonicalMatchHighDL, Drivers.canonicalMutationEntryForGene(drivers, "gene 1"));
        assertNull(Drivers.canonicalMutationEntryForGene(drivers, "gene 2"));
    }

}