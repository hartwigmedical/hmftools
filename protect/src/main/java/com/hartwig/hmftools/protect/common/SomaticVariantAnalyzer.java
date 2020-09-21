package com.hartwig.hmftools.protect.common;

import java.io.IOException;
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
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SomaticVariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private SomaticVariantAnalyzer() {

    }

    @NotNull
    public static SomaticVariantAnalysis analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @NotNull List<DriverCatalog> driverCatalog) throws IOException {

        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        List<SomaticVariant> variantsToReport = variants.stream().filter(x -> x.reported()).collect(Collectors.toList());

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, driverCatalog);

    }
}
