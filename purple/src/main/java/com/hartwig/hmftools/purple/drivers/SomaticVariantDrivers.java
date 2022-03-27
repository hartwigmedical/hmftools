package com.hartwig.hmftools.purple.drivers;

import java.util.List;
import java.util.Map;
import java.util.Optional;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverImpact;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.ReportablePredicate;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class SomaticVariantDrivers
{
    private final DriverGenePanel mGenePanel;

    private final List<SomaticVariant> mTsgVariants;
    private final List<SomaticVariant> mOncoVariants;
    private final Map<VariantType,Integer> mVariantTypeCounts;
    private final Map<VariantType,Integer> mVariantTypeCountsBiallelic;

    private final OncoDrivers mOncoDrivers;
    private final TsgDrivers mTsgDrivers;

    public SomaticVariantDrivers(final DriverGenePanel panel)
    {
        mGenePanel = panel;

        mTsgVariants = Lists.newArrayList();
        mOncoVariants = Lists.newArrayList();
        mVariantTypeCounts = Maps.newHashMap();
        mVariantTypeCountsBiallelic = Maps.newHashMap();

        mOncoDrivers = new OncoDrivers(panel);
        mTsgDrivers = new TsgDrivers(panel);
    }

    public boolean checkSomaticVariant(final SomaticVariant variant)
    {
        // return true if the variant is a reportable TSG or onocogene
        if(variant.isFiltered())
            return false;

        mVariantTypeCounts.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);

        if(variant.biallelic())
        {
            mVariantTypeCountsBiallelic.compute(variant.type(), (key, oldValue) -> Optional.ofNullable(oldValue).orElse(0) + 1);
        }

        if(mOncoDrivers.checkVariant(variant))
            return true;

        if(mTsgDrivers.checkVariant(variant))
            return true;

        return false;
    }

    protected static boolean isReportable(final ReportablePredicate predicate, final SomaticVariant variant)
    {
        final VariantImpact variantImpact = variant.variantImpact();

        return predicate.test(
                variantImpact.CanonicalGeneName, variant.type(), variant.decorator().repeatCount(), variant.isHotspot(),
                variantImpact.CanonicalCodingEffect, variantImpact.CanonicalEffect);
    }

    public List<DriverCatalog> buildCatalog(final Map<String,List<GeneCopyNumber>> geneCopyNumberMap)
    {
        final List<DriverCatalog> result = Lists.newArrayList();

        result.addAll(mOncoDrivers.findDrivers(geneCopyNumberMap, mVariantTypeCounts));
        result.addAll(mTsgDrivers.findDrivers(geneCopyNumberMap, mVariantTypeCounts, mVariantTypeCountsBiallelic));

        return result;
    }

    protected static Map<DriverImpact,Integer> groupByImpact(final List<SomaticVariant> variants)
    {
        Map<DriverImpact,Integer> map = Maps.newHashMap();

        for(SomaticVariant variant : variants)
        {
            DriverImpact driverImpact = DriverImpact.select(variant.type(), variant.variantImpact().CanonicalCodingEffect);
            Integer count = map.get(driverImpact);
            map.put(driverImpact, count != null ? count + 1 : 1);
        }

        return map;
    }
}
