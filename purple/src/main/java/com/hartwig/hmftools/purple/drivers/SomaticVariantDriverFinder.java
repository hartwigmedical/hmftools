package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.isReportable;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.purple.DriverGeneResource;
import com.hartwig.hmftools.common.driver.panel.ReportablePredicate;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.DriverSourceData;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public abstract class SomaticVariantDriverFinder
{
    protected final DriverCategory mCategory;
    protected final ReportablePredicate mReportablePredicate;
    protected final Map<String,DndsDriverGeneLikelihood> mLikelihoodsByGene;
    protected final List<SomaticVariant> mReportableVariants;

    public static final DndsDriverGeneLikelihood NO_GENE_DNDS_LIKELIHOOD = ImmutableDndsDriverGeneLikelihood.builder()
        .gene("")
        .indel(ImmutableDndsDriverImpactLikelihood.builder().driversPerSample(0).passengersPerMutation(0).build())
        .missense(ImmutableDndsDriverImpactLikelihood.builder().driversPerSample(0).passengersPerMutation(0).build())
        .splice(ImmutableDndsDriverImpactLikelihood.builder().driversPerSample(0).passengersPerMutation(0).build())
        .nonsense(ImmutableDndsDriverImpactLikelihood.builder().driversPerSample(0).passengersPerMutation(0).build())
        .build();

    public SomaticVariantDriverFinder(final DriverGeneResource genePanel, final DriverCategory category)
    {
        mCategory = category;
        mLikelihoodsByGene = category == DriverCategory.ONCO ? genePanel.OncoLikelihoodMap : genePanel.TsgLikelihoodMap;
        mReportablePredicate = new ReportablePredicate(category, genePanel.DriverGeneList);
        mReportableVariants = Lists.newArrayList();
    }
    
    public boolean checkVariant(final SomaticVariant variant)
    {
        if(isReportable(mReportablePredicate, variant))
        {
            mReportableVariants.add(variant);
            return true;
        }
        
        return false;
    }

    public boolean hasVariant(final SomaticVariant variant) { return mReportableVariants.contains(variant); }
    public void addVariant(final SomaticVariant variant) { mReportableVariants.add(variant); }

    public abstract List<DriverCatalog> findDrivers(
            final Map<String,List<GeneCopyNumber>> geneCopyNumberMap, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> variantTypeCountsBiallelic, final List<DriverSourceData> driverSourceData);
}
