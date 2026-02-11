package com.hartwig.hmftools.purple.drivers;

import static com.hartwig.hmftools.common.driver.LikelihoodMethod.DNDS;
import static com.hartwig.hmftools.common.driver.LikelihoodMethod.SPLICE_REGION;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.hasTranscriptCodingEffect;
import static com.hartwig.hmftools.purple.drivers.SomaticVariantDrivers.isCandidateReportable;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.LikelihoodMethod;
import com.hartwig.hmftools.common.driver.dnds.DndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverGeneLikelihood;
import com.hartwig.hmftools.common.driver.dnds.ImmutableDndsDriverImpactLikelihood;
import com.hartwig.hmftools.common.purple.ReportedStatus;
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
    protected final List<SomaticVariant> mDriverVariants;

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
        mDriverVariants = Lists.newArrayList();
    }
    
    public boolean checkVariant(final SomaticVariant variant)
    {
        if(isCandidateReportable(mReportablePredicate, variant))
        {
            mDriverVariants.add(variant);

            return mReportablePredicate.isReportable(variant.variantImpact(), variant.isHotspot());
        }
        
        return false;
    }

    public boolean hasVariant(final SomaticVariant variant) { return mDriverVariants.contains(variant); }
    public void addVariant(final SomaticVariant variant) { mDriverVariants.add(variant); }

    protected List<DriverCatalog> findDrivers(
            final Map<String,List<GeneCopyNumber>> geneCopyNumberMap, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> variantTypeCountsBiallelic, final List<DriverSourceData> driverSourceData)
    {
        List<DriverCatalog> driverCatalog = Lists.newArrayList();

        Map<String,List<SomaticVariant>> codingVariants = mDriverVariants.stream().collect(Collectors.groupingBy(SomaticVariant::gene));

        for(Map.Entry<String,List<SomaticVariant>> entry : codingVariants.entrySet())
        {
            String gene = entry.getKey();
            List<SomaticVariant> geneVariants = entry.getValue();

            DndsDriverGeneLikelihood dndsLikelihood = mLikelihoodsByGene.containsKey(gene) ?
                    mLikelihoodsByGene.get(gene) : NO_GENE_DNDS_LIKELIHOOD;

            List<GeneCopyNumber> geneCopyNumbers = geneCopyNumberMap.get(gene);

            if(geneCopyNumbers == null)
                continue;

            if(!mReportablePredicate.isDriverGene(geneCopyNumbers.get(0).GeneName))
                continue;

            for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
            {
                if(geneCopyNumbers.size() == 1
                || geneVariants.stream().anyMatch(x -> hasTranscriptCodingEffect(x.variantImpact(), x.type(), geneCopyNumber.TransName)))
                {
                    // check the driver-gene-panel reportable flag based on the type of mutation
                    boolean hasCodingImpact = geneVariants.stream()
                            .anyMatch(x -> x.isHotspot() || hasTranscriptCodingEffect(x.variantImpact(), x.type(), geneCopyNumber.TransName));

                    boolean isReportable = hasCodingImpact
                            && geneVariants.stream().anyMatch(x -> mReportablePredicate.isReportable(x.variantImpact(), x.isHotspot()));

                    ReportedStatus reportedStatus = isReportable ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;

                    LikelihoodMethod likelihoodMethod = hasCodingImpact ? DNDS : SPLICE_REGION;

                    DriverCatalog driverRecord = createDriverCatalog(
                            geneVariants, variantTypeCounts, variantTypeCountsBiallelic, geneCopyNumber,
                            dndsLikelihood, likelihoodMethod, reportedStatus);

                    driverCatalog.add(driverRecord);

                    driverSourceData.add(new DriverSourceData(driverRecord, geneVariants.get(0)));

                }
            }
        }

        return driverCatalog;
    }

    public abstract DriverCatalog createDriverCatalog(
            final List<SomaticVariant> geneVariants, final Map<VariantType,Integer> variantTypeCounts,
            final Map<VariantType,Integer> biallelicCounts, final GeneCopyNumber geneCopyNumber,
            final DndsDriverGeneLikelihood likelihood, final LikelihoodMethod likelihoodMethod, final ReportedStatus reportedStatus);
}
