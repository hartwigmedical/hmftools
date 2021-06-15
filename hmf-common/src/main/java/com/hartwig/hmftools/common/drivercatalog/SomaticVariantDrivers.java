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

public class SomaticVariantDrivers
{
    private final DriverGenePanel mGenePanel;

    private final List<SomaticVariant> mTsgVariants;
    private final List<SomaticVariant> mOncoVariants;
    private final Map<VariantType, Long> mVariantTypeCounts;
    private final Map<VariantType, Long> mVariantTypeCountsBiallelic;

    private final Predicate<SomaticVariant> mOncoPredicate;
    private final Predicate<SomaticVariant> mTsgPredicate;

    public SomaticVariantDrivers(@NotNull final DriverGenePanel panel)
    {
        mGenePanel = panel;

        mTsgVariants = Lists.newArrayList();
        mOncoVariants = Lists.newArrayList();
        mVariantTypeCounts = Maps.newHashMap();
        mVariantTypeCountsBiallelic = Maps.newHashMap();

        mOncoPredicate = new ReportablePredicate(DriverCategory.ONCO, panel);
        mTsgPredicate = new ReportablePredicate(DriverCategory.TSG, panel);
    }

    public boolean add(@NotNull final SomaticVariant variant)
    {
        if(!variant.isFiltered())
        {
            mVariantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            if(variant.biallelic())
            {
                mVariantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0L) + 1);
            }

            if(mOncoPredicate.test(variant))
            {
                mOncoVariants.add(variant);
                return true;
            }

            if(mTsgPredicate.test(variant))
            {
                mTsgVariants.add(variant);
                return true;
            }
        }

        return false;
    }

    @NotNull
    public List<DriverCatalog> build(@NotNull final List<GeneCopyNumber> geneCopyNumbers)
    {
        final OncoDrivers oncoDrivers = new OncoDrivers(mGenePanel);
        final TsgDrivers tsgDrivers = new TsgDrivers(mGenePanel);

        final List<DriverCatalog> result = Lists.newArrayList();
        result.addAll(oncoDrivers.drivers(mOncoVariants, geneCopyNumbers, mVariantTypeCounts));
        result.addAll(tsgDrivers.drivers(mTsgVariants, geneCopyNumbers, mVariantTypeCounts, mVariantTypeCountsBiallelic));

        return result;
    }
}
