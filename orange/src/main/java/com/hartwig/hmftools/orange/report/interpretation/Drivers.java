package com.hartwig.hmftools.orange.report.interpretation;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class Drivers {

    private static final Set<DriverType> VARIANT_DRIVER_TYPES = Sets.newHashSet(DriverType.MUTATION, DriverType.GERMLINE_MUTATION);

    private Drivers() {
    }

    @Nullable
    public static DriverCatalog variantEntryForGene(@NotNull List<DriverCatalog> drivers, @NotNull String geneToFind) {
        DriverCatalog highest = null;
        for (DriverCatalog driver : drivers) {
            if (VARIANT_DRIVER_TYPES.contains(driver.driver()) && driver.gene().equals(geneToFind)) {
                if (highest == null || driver.driverLikelihood() > highest.driverLikelihood()) {
                    highest = driver;
                }
            }
        }

        return highest;
    }
}
