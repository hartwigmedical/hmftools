package com.hartwig.hmftools.purple.fitting;

import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.purple.PurpleUtils.PPL_LOGGER;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_TOTAL_SOMATIC_VAR_ALLELE_READ_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.MIN_TOTAL_SV_FRAGMENT_COUNT;
import static com.hartwig.hmftools.purple.PurpleConstants.NO_TUMOR_BAF_TOTAL;
import static com.hartwig.hmftools.purple.PurpleConstants.NO_TUMOR_DEPTH_RATIO_MAX;
import static com.hartwig.hmftools.purple.PurpleConstants.NO_TUMOR_DEPTH_RATIO_MIN;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.FittedPurity;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.purple.PurpleConfig;
import com.hartwig.hmftools.purple.ReferenceData;
import com.hartwig.hmftools.purple.SampleData;
import com.hartwig.hmftools.purple.SomaticFitConfig;
import com.hartwig.hmftools.purple.fittingsnv.SomaticPurityFitter;
import com.hartwig.hmftools.purple.region.ObservedRegion;
import com.hartwig.hmftools.purple.somatic.SomaticVariant;

public class VariantPurityFitter
{
    private final ReferenceData mReferenceData;
    private final SampleData mSampleData;
    private final SomaticPurityFitter mSomaticPurityFitter;

    // interim state
    private List<SomaticVariant> mFittingSomatics;

    private int mSvHotspotCount;
    private int mSvFragmentReadCount;
    private int mSomaticHotspotCount;
    private int mAlleleReadCountTotal;

    private boolean mHasTumor;

    public VariantPurityFitter(final PurpleConfig config, final ReferenceData referenceData, final SampleData sampleData)
    {
        mReferenceData = referenceData;
        mSampleData = sampleData;

        mSvHotspotCount = 0;
        mSvFragmentReadCount = 0;
        mSomaticHotspotCount = 0;
        mAlleleReadCountTotal = 0;
        mHasTumor = false;

        mFittingSomatics = Lists.newArrayList();

        mSomaticPurityFitter = new SomaticPurityFitter(
                config.SomaticFitting.MinPeakVariants, config.SomaticFitting.MinTotalVariants,
                sampleData.Amber.minSomaticTotalReadCount(), sampleData.Amber.maxSomaticTotalReadCount(),
                config.Fitting.MinPurity, config.Fitting.MaxPurity);
    }

    public boolean hasTumor() { return mHasTumor; }
    public List<SomaticVariant> fittingSomatics() { return mFittingSomatics; }

    public void setState(final List<ObservedRegion> observedRegions)
    {
        mFittingSomatics.addAll(SomaticPurityFitter.findFittingVariants(mSampleData.SomaticCache.variants(), observedRegions));

        if(!mFittingSomatics.isEmpty())
        {
            PPL_LOGGER.debug("somatic fitting variants({})", mFittingSomatics.size());
        }

        setSvSummary(mSampleData.SvCache.variants());
        setSomaticSummary(mFittingSomatics);

        if(mSomaticHotspotCount > 0 || mAlleleReadCountTotal >= MIN_TOTAL_SOMATIC_VAR_ALLELE_READ_COUNT)
        {
            PPL_LOGGER.debug("tumor evidence: somaticHotspotCount({}) alleleReadCountTotal({})",
                    mSomaticHotspotCount, mAlleleReadCountTotal);
            mHasTumor = true;
            return;
        }

        if(mSvHotspotCount > 0 || mSvFragmentReadCount >= MIN_TOTAL_SV_FRAGMENT_COUNT)
        {
            PPL_LOGGER.debug("tumor evidence: svHotspotCount({}) svFragmentReadCount({})",
                    mSvHotspotCount, mSvFragmentReadCount);
            mHasTumor = true;
            return;
        }

        int tumorEvidenceBafCountTotal = observedRegions.stream()
                .filter(x -> x.germlineStatus() == DIPLOID)
                .filter(x -> x.observedTumorRatio() < NO_TUMOR_DEPTH_RATIO_MIN || x.observedTumorRatio() > NO_TUMOR_DEPTH_RATIO_MAX)
                .mapToInt(x -> x.bafCount())
                .sum();

        if(tumorEvidenceBafCountTotal >= NO_TUMOR_BAF_TOTAL)
        {
            PPL_LOGGER.debug("tumor evidence: tumorEvidenceBafCountTotal({})", tumorEvidenceBafCountTotal);
            mHasTumor = true;
            return;
        }

        mHasTumor = false;
    }

    public FittedPurity calcSomaticFit(final List<FittedPurity> diploidCandidates, final List<PurpleCopyNumber> copyNumbers)
    {
        return mSomaticPurityFitter.fromSomatics(mFittingSomatics, diploidCandidates, copyNumbers);
    }

    public FittedPurity tumorOnlySomaticFit(final List<FittedPurity> allCandidates)
    {
        return mSomaticPurityFitter.fromTumorOnlySomatics(mReferenceData.DriverGenes, mFittingSomatics, allCandidates);
    }

    public static boolean somaticFitIsWorse(final FittedPurity lowestScore, final FittedPurity somaticFit, final SomaticFitConfig config)
    {
        double lowestPurity = lowestScore.purity();
        double somaticPurity = somaticFit.purity();

        return Doubles.lessThan(lowestPurity, config.MinSomaticPurity)
            && Doubles.lessThan(somaticPurity, config.MinSomaticPurity)
            && Doubles.greaterThan( somaticPurity, lowestPurity);
    }

    private void setSvSummary(final List<StructuralVariant> variants)
    {
        for(StructuralVariant variant : variants)
        {
            if(variant.isFiltered())
                continue;

            if(variant.hotspot())
                mSvHotspotCount++;

            Integer startTumorVariantFragmentCount = variant.start().tumorVariantFragmentCount();
            if(variant.end() != null && startTumorVariantFragmentCount != null)
            {
                mSvFragmentReadCount += startTumorVariantFragmentCount;
            }
        }
    }

    private void setSomaticSummary(final List<SomaticVariant> somatics)
    {
        for(SomaticVariant variant : somatics)
        {
            if(!variant.isPass())
                continue;

            if(variant.isHotspot())
                mSomaticHotspotCount++;

            if(variant.type() == VariantType.SNP) // they all should be
            {
                mAlleleReadCountTotal += variant.alleleReadCount();
            }
        }
    }
}
