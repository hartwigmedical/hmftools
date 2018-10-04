package com.hartwig.hmftools.patientreporter.algo;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.dnds.DndsDriverGeneLikelihoodSupplier;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DriverProbabilityModel {

    @NotNull
    private final List<DriverCatalog> driverCatalogList;

    @NotNull
    public static List<DriverCatalog> createDriverCatalogForSomaticVariants(@NotNull List<EnrichedSomaticVariant> variants) {
        final List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(DndsDriverGeneLikelihoodSupplier.oncoLikelihood(), variants));
        driverCatalog.addAll(TsgDrivers.drivers(DndsDriverGeneLikelihoodSupplier.tsgLikelihood(), variants));

        return driverCatalog;
    }

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
