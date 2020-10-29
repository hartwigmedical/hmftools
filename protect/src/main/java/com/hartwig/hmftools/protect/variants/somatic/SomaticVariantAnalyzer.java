package com.hartwig.hmftools.protect.variants.somatic;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.jetbrains.annotations.NotNull;

@Deprecated
public final class SomaticVariantAnalyzer {

    private SomaticVariantAnalyzer() {
    }

    @NotNull
    public static List<DriverSomaticVariant> run(@NotNull List<SomaticVariant> variants, @NotNull List<DriverCatalog> purpleDriverCatalog) {
        List<DriverSomaticVariant> driverSomaticVariants = Lists.newArrayList();

        for (SomaticVariant variant : variants) {
            if (variant.reported()) {
                DriverCatalog entry = catalogEntryForGene(purpleDriverCatalog, variant.gene());
                driverSomaticVariants.add(ImmutableDriverSomaticVariant.builder()
                        .variant(variant)
                        .driverLikelihood(entry.driverLikelihood())
                        .build());
            }
        }

        return driverSomaticVariants;
    }

    @NotNull
    private static DriverCatalog catalogEntryForGene(@NotNull List<DriverCatalog> driverCatalogList, @NotNull String gene) {
        for (DriverCatalog entry : driverCatalogList) {
            if (entry.gene().equals(gene)) {
                return entry;
            }
        }

        throw new IllegalStateException("Could not find driver catalog entry for " + gene);
    }
}
