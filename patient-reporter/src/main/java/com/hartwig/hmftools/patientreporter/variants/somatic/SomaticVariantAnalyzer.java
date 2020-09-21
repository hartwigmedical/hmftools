package com.hartwig.hmftools.patientreporter.variants.somatic;

import java.util.List;

import java.util.stream.Collectors;

import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

public final class SomaticVariantAnalyzer {

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<SomaticVariant> variants, @NotNull List<DriverCatalog> driverCatalog) {
        List<SomaticVariant> variantsToReport = variants.stream().filter(x -> x.reported()).collect(Collectors.toList());

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, driverCatalog);
    }

}
