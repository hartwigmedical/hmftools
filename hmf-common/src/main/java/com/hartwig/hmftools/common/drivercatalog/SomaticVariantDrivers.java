package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class SomaticVariantDrivers {

    private final DriverGenePanel genePanel;

    private final List<SomaticVariant> tsgVariants = Lists.newArrayList();
    private final List<SomaticVariant> oncoVariants = Lists.newArrayList();
    private final Map<VariantType, Long> variantTypeCounts = Maps.newHashMap();
    private final Map<VariantType, Long> variantTypeCountsBiallelic = Maps.newHashMap();

    private final Predicate<SomaticVariant> oncoPredicate;
    private final Predicate<SomaticVariant> tsgPredicate;

    public SomaticVariantDrivers(@NotNull final DriverGenePanel panel) {
        this.genePanel = panel;
        oncoPredicate = new ReportablePredicate(DriverCategory.ONCO, panel);
        tsgPredicate = new ReportablePredicate(DriverCategory.TSG, panel);
    }

    public void add(@NotNull final SomaticVariant variant) {
        if (!variant.isFiltered()) {
            variantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            if (variant.biallelic()) {
                variantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            }

            if (oncoPredicate.test(variant)) {
                oncoVariants.add(variant);
            }

            if (tsgPredicate.test(variant)) {
                tsgVariants.add(variant);
            }
        }
    }

    @NotNull
    public List<DriverCatalog> build(@NotNull final List<GeneCopyNumber> geneCopyNumbers) {
        final OncoDrivers oncoDrivers = new OncoDrivers(genePanel);
        final TsgDrivers tsgDrivers = new TsgDrivers(genePanel);

        final List<DriverCatalog> result = Lists.newArrayList();
        result.addAll(oncoDrivers.drivers(oncoVariants, geneCopyNumbers, variantTypeCounts));
        result.addAll(tsgDrivers.drivers(tsgVariants, geneCopyNumbers, variantTypeCounts, variantTypeCountsBiallelic));

        return result;
    }

}
