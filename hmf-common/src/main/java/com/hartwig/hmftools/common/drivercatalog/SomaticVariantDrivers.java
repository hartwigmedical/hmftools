package com.hartwig.hmftools.common.drivercatalog;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class SomaticVariantDrivers {

    private final List<SomaticVariant> tsgVariants = Lists.newArrayList();
    private final List<SomaticVariant> oncoVariants = Lists.newArrayList();
    private final Map<VariantType, Long> variantTypeCounts = Maps.newHashMap();
    private final Map<VariantType, Long> variantTypeCountsBiallelic = Maps.newHashMap();
    private final Map<VariantType, Long> variantTypeCountsNonBiallelic = Maps.newHashMap();
    private final Map<String, DndsDriverGeneLikelihood> tsgLikelihood;
    private final Map<String, DndsDriverImpactLikelihood> oncoLikelihood;

    private final Predicate<SomaticVariant> oncoPredicate;
    private final Predicate<SomaticVariant> tsgPredicate;

    public SomaticVariantDrivers(@NotNull final DriverGenePanel panel) {
        tsgLikelihood = panel.tsgLikelihood();
        oncoLikelihood = panel.oncoLikelihood();
        oncoPredicate = OncoDrivers.oncoVariant(panel.oncoGenes());
        tsgPredicate = TsgDrivers.tsgVariant(panel.tsGenes());
    }

    public void add(@NotNull final SomaticVariant variant) {
        if (!variant.isFiltered()) {
            variantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            if (variant.biallelic()) {
                variantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            } else {
                variantTypeCountsNonBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
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

        final List<DriverCatalog> result = Lists.newArrayList();
        result.addAll(OncoDrivers.drivers(oncoLikelihood, oncoVariants, geneCopyNumbers, variantTypeCounts));
        result.addAll(TsgDrivers.drivers(tsgLikelihood,
                tsgVariants,
                geneCopyNumbers,
                variantTypeCounts,
                variantTypeCountsBiallelic,
                variantTypeCountsNonBiallelic));

        return result;
    }

}
