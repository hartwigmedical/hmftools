package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverProbabilityModel {

    @NotNull
    private final List<DriverCatalog> driverCatalogList;

    public DriverProbabilityModel(@NotNull List<DriverCatalog> driverCatalogList) {
        this.driverCatalogList = driverCatalogList;

    }

    @Nullable
    public DriverCatalog catalogForVariant(@NotNull SomaticVariant variant) {
        // TODO (KODU): Filter on snv-like catalog entries only.
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(variant.gene())) {
                return entry;
            }
        }
        return null;
    }
}
