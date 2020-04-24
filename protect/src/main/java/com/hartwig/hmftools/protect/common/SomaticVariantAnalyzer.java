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
import com.hartwig.hmftools.common.drivercatalog.OncoDrivers;
import com.hartwig.hmftools.common.drivercatalog.TsgDrivers;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteIndels;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalLoad;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class SomaticVariantAnalyzer {

    private static final Logger LOGGER = LogManager.getLogger(SomaticVariantAnalyzer.class);

    private static final List<CodingEffect> TSG_CODING_EFFECTS_TO_REPORT =
            Lists.newArrayList(CodingEffect.NONSENSE_OR_FRAMESHIFT, CodingEffect.MISSENSE, CodingEffect.SPLICE);

    private static final Set<String> ONCO_GENES_WITH_SPLICE_EFFECTS = Sets.newHashSet("MET", "JAK2");

    private static final List<CodingEffect> ONCO_CODING_EFFECTS_TO_REPORT = Lists.newArrayList(CodingEffect.MISSENSE);

    private SomaticVariantAnalyzer() {

    }

    @NotNull
    public static SomaticVariantAnalysis analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf,
            @NotNull List<GeneCopyNumber> exomeGeneCopyNumbers, @NotNull DriverGeneView driverGeneView) throws IOException {

        List<SomaticVariant> variants = SomaticVariantFactory.passOnlyInstance().fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info("Loaded {} PASS somatic variants from {}", variants.size(), somaticVariantVcf);

        List<SomaticVariant> variantsToReport = variants.stream().filter(includeFilter(driverGeneView)).collect(Collectors.toList());

        List<DriverCatalog> driverCatalog = Lists.newArrayList();
        driverCatalog.addAll(OncoDrivers.drivers(variants, exomeGeneCopyNumbers));
        driverCatalog.addAll(TsgDrivers.drivers(variants, exomeGeneCopyNumbers));

        // Check that we miss no drivers
        for (DriverCatalog driver : driverCatalog) {
            boolean reported = false;
            for (SomaticVariant variant : variantsToReport) {
                if (variant.gene().equals(driver.gene())) {
                    reported = true;
                }
            }
            if (!reported) {
                LOGGER.warn("Driver {} not added to reported somatic variants", driver.gene());
            }
        }

        return ImmutableSomaticVariantAnalysis.of(variantsToReport,
                driverCatalog);

    }

    @NotNull
    private static Predicate<SomaticVariant> includeFilter(@NotNull DriverGeneView driverGeneView) {
        return variant -> {
            if (variant.isFiltered()) {
                return false;
            }

            if (driverGeneView.oncoDriverGenes().contains(variant.gene()) || driverGeneView.tsgDriverGenes().contains(variant.gene())) {
                if (variant.isHotspot()) {
                    return true;
                }

                DriverCategory category = driverGeneView.category(variant.gene());
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
