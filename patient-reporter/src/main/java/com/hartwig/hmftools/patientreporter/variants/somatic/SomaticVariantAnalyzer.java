package com.hartwig.hmftools.patientreporter.variants.somatic;

import java.util.List;
import java.util.Set;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.SomaticVariantDrivers;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class SomaticVariantAnalyzer {

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<SomaticVariant> variants, @NotNull DriverGenePanel driverGenePanel,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers) {
        List<SomaticVariant> variantsToReport = variants.stream().filter(x -> x.reported()).collect(Collectors.toList());

        SomaticVariantDrivers drivers = new SomaticVariantDrivers(driverGenePanel);
        variants.forEach(drivers::add);
        List<DriverCatalog> driverCatalog = drivers.build(exomeGeneCopyNumbers);

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, driverCatalog);
    }

}
