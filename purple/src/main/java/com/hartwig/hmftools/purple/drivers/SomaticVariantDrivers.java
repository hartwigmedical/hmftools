package com.hartwig.hmftools.purple.drivers;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.function.Predicate;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;

import org.jetbrains.annotations.NotNull;

public class SomaticVariantDrivers
{
    private final DriverGenePanel mGenePanel;

    private final List<SomaticVariant> mTsgVariants;
    private final List<SomaticVariant> mOncoVariants;
    private final Map<VariantType,Integer> mVariantTypeCounts;
    private final Map<VariantType,Integer> mVariantTypeCountsBiallelic;

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
            mVariantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);
            if(variant.biallelic())
            {
                mVariantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);
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

    public List<DriverCatalog> build(final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        final OncoDrivers oncoDrivers = new OncoDrivers(mGenePanel);
        final TsgDrivers tsgDrivers = new TsgDrivers(mGenePanel);

        final List<DriverCatalog> result = Lists.newArrayList();

        result.addAll(oncoDrivers.drivers(mOncoVariants, geneCopyNumberMap, mVariantTypeCounts));
        result.addAll(tsgDrivers.drivers(mTsgVariants, geneCopyNumberMap, mVariantTypeCounts, mVariantTypeCountsBiallelic));

        return result;
    }
}
