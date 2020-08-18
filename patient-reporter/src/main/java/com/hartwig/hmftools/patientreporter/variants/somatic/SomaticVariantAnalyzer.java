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

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final Set<String> ONCO_GENES_WITH_SPLICE_EFFECTS = Sets.newHashSet("MET", "JAK2");

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static SomaticVariantAnalysis run(@NotNull List<SomaticVariant> variants, @NotNull DriverGenePanel driverGenePanel,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers) {
        List<SomaticVariant> variantsToReport = variants.stream().filter(includeFilter(driverGenePanel)).collect(Collectors.toList());

        SomaticVariantDrivers drivers = new SomaticVariantDrivers(driverGenePanel);
        variants.forEach(drivers::add);
        List<DriverCatalog> driverCatalog = drivers.build(exomeGeneCopyNumbers);

        // Check that we miss no driver genes.
        for (DriverCatalog driver : driverCatalog) {
            boolean reported = false;
            for (SomaticVariant variant : variantsToReport) {
                if (variant.gene().equals(driver.gene())) {
                    reported = true;
                }
            }
            if (!reported) {
                LOGGER.warn("Driver gene '{}' not added to reported somatic variants!", driver.gene());
            }
        }

        return ImmutableSomaticVariantAnalysis.of(variantsToReport, driverCatalog);
    }

    @NotNull
    private static Predicate<SomaticVariant> includeFilter(@NotNull DriverGenePanel driverGenePanel) {
        return variant -> {
            if (variant.isFiltered()) {
                return false;
            }

            if (driverGenePanel.oncoGenes().contains(variant.gene()) || driverGenePanel.tsGenes().contains(variant.gene())) {
                if (variant.isHotspot()) {
                    return true;
                }

                DriverCategory category = driverGenePanel.category(variant.gene());
                if (category == null) {
                    throw new IllegalStateException("Driver category not known for driver gene: " + variant.gene());
                }

                CodingEffect effect = variant.canonicalCodingEffect();
                if (category == DriverCategory.TSG) {
                    return TSG_CODING_EFFECTS_TO_REPORT.contains(effect);
                } else {
                    assert category == DriverCategory.ONCO;
                    if (ONCO_GENES_WITH_SPLICE_EFFECTS.contains(variant.gene())) {
                        return ONCO_CODING_EFFECTS_TO_REPORT.contains(effect) || effect == CodingEffect.SPLICE;
                    } else {
                        return ONCO_CODING_EFFECTS_TO_REPORT.contains(effect);
                    }
                }
            }
            return false;
        };
    }
}
