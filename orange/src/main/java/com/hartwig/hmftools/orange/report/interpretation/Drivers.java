package com.hartwig.hmftools.orange.report.interpretation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Drivers {

    private static final Set<DriverType> MUTATION_DRIVER_TYPES = Sets.newHashSet(DriverType.MUTATION, DriverType.GERMLINE_MUTATION);

    private Drivers() {
    }

    @NotNull
    public static List<DriverCatalog> nonCanonicalMutationEntries(@NotNull List<DriverCatalog> drivers) {
        List<DriverCatalog> nonCanonicalVariantEntries = Lists.newArrayList();
        for (DriverCatalog driver : drivers) {
            if (MUTATION_DRIVER_TYPES.contains(driver.driver()) && !driver.isCanonical()) {
                nonCanonicalVariantEntries.add(driver);
            }
        }
        return nonCanonicalVariantEntries;
    }

    @Nullable
    public static DriverCatalog canonicalMutationEntryForGene(@NotNull List<DriverCatalog> drivers, @NotNull String geneToFind) {
        DriverCatalog highest = null;
        for (DriverCatalog driver : drivers) {
            if (MUTATION_DRIVER_TYPES.contains(driver.driver()) && driver.gene().equals(geneToFind) && driver.isCanonical()) {
                if (highest == null || driver.driverLikelihood() > highest.driverLikelihood()) {
                    highest = driver;
                }
            }
        }

        return highest;
    }
}
