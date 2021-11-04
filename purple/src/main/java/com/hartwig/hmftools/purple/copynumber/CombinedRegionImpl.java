package com.hartwig.hmftools.purple.copynumber;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.purple.copynumber.CopyNumberMethod;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.purple.region.GermlineStatus;
import com.hartwig.hmftools.common.purple.region.ModifiableFittedRegion;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;
import com.hartwig.hmftools.common.utils.Doubles;

import org.jetbrains.annotations.NotNull;

class CombinedRegionImpl implements CombinedRegion
{
    private final ModifiableFittedRegion mCombined;

    private CopyNumberMethod mCopyNumberMethod = CopyNumberMethod.UNKNOWN;
    private boolean mInferredBAF;
    private final List<FittedRegion> mRegions = Lists.newArrayList();
    private int mUnweightedCount = 1;
    private final boolean mIsBafWeighted;
    private SegmentSupport mGermlineEndSupport = SegmentSupport.UNKNOWN;

    CombinedRegionImpl(@NotNull final FittedRegion region)
    {
        this(true, region);
    }

    CombinedRegionImpl(boolean isBafWeighed, @NotNull final FittedRegion region)
    {
        mCombined = ModifiableFittedRegion.create().from(region);
        mIsBafWeighted = isBafWeighed;

        if(region.germlineStatus() != GermlineStatus.DIPLOID)
        {
            clearBAFValues();
        }
        mRegions.add(region);
    }

    @NotNull
    @Override
    public String chromosome()
    {
        return mCombined.chromosome();
    }

    @Override
    public long start()
    {
        return mCombined.start();
    }

    @Override
    public long end()
    {
        return mCombined.end();
    }

    public boolean isInferredBAF()
    {
        return mInferredBAF;
    }

    @NotNull
    public List<FittedRegion> regions()
    {
        return mRegions;
    }

    public double tumorCopyNumber()
    {
        return mCombined.tumorCopyNumber();
    }

    public double tumorBAF()
    {
        return mCombined.tumorBAF();
    }

    public int bafCount()
    {
        return region().bafCount();
    }

    @NotNull
    public CopyNumberMethod copyNumberMethod()
    {
        return mCopyNumberMethod;
    }

    public boolean isProcessed()
    {
        return mCopyNumberMethod != CopyNumberMethod.UNKNOWN;
    }

    public void setCopyNumberMethod(@NotNull CopyNumberMethod copyNumberMethod)
    {
        mCopyNumberMethod = copyNumberMethod;
    }

    @NotNull
    public FittedRegion region()
    {
        return mCombined;
    }

    @NotNull
    @Override
    public SegmentSupport germlineEndSupport()
    {
        return mGermlineEndSupport;
    }

    public void setGermlineEndSupport(@NotNull final SegmentSupport germlineEndSupport)
    {
        mGermlineEndSupport = germlineEndSupport;
    }

    public void extend(@NotNull final FittedRegion region)
    {
        mCombined.setStart(Math.min(mCombined.start(), region.start()));
        mCombined.setEnd(Math.max(mCombined.end(), region.end()));
        mCombined.setMinStart(Math.min(region.minStart(), region().minStart()));
        mCombined.setMaxStart(Math.min(region.maxStart(), region().maxStart()));

        if(region.start() <= mCombined.start())
        {
            mRegions.add(0, region);
            mCombined.setSupport(region.support());
            mCombined.setRatioSupport(region.ratioSupport());
        }
        else
        {
            mRegions.add(region);
        }
    }

    public void extendWithUnweightedAverage(@NotNull final FittedRegion region)
    {
        extend(region);
        applyDepthWindowCountWeights(region, mUnweightedCount, 1);
        applyBafCountWeights(region, mUnweightedCount, 1);
        mUnweightedCount++;
    }

    public void extendWithWeightedAverage(@NotNull final FittedRegion region)
    {
        mCombined.setGermlineStatus(GermlineStatus.DIPLOID);

        int currentWeight = region().depthWindowCount();
        int newWeight = region.depthWindowCount();

        applyDepthWindowCountWeights(region, currentWeight, newWeight);

        if(mIsBafWeighted && (mCombined.bafCount() > 0 || region.bafCount() > 0))
        {
            currentWeight = mCombined.bafCount();
            newWeight = region.bafCount();
        }
        applyBafCountWeights(region, currentWeight, newWeight);

        extend(region);
    }

    private void applyBafCountWeights(final FittedRegion region, long currentWeight, long newWeight)
    {
        if(!Doubles.isZero(region.observedBAF()))
        {
            mCombined.setObservedBAF(weightedAverage(currentWeight, mCombined.observedBAF(), newWeight, region.observedBAF()));
        }

        if(!Doubles.isZero(region.tumorBAF()))
        {
            mCombined.setTumorBAF(weightedAverage(currentWeight, mCombined.tumorBAF(), newWeight, region.tumorBAF()));
        }

        mCombined.setBafCount(mCombined.bafCount() + region.bafCount());
    }

    private void applyDepthWindowCountWeights(final FittedRegion region, long currentWeight, long newWeight)
    {

        if(!Doubles.isZero(region.tumorCopyNumber()))
        {
            mCombined.setTumorCopyNumber(weightedAverage(currentWeight, mCombined.tumorCopyNumber(), newWeight, region.tumorCopyNumber()));
        }

        if(!Doubles.isZero(region.refNormalisedCopyNumber()))
        {
            mCombined.setRefNormalisedCopyNumber(weightedAverage(currentWeight,
                    mCombined.refNormalisedCopyNumber(),
                    newWeight,
                    region.refNormalisedCopyNumber()));
        }

        if(!Doubles.isZero(region.gcContent()))
        {
            mCombined.setGcContent(weightedAverage(currentWeight, mCombined.gcContent(), newWeight, region.gcContent()));
        }

        mCombined.setDepthWindowCount(mCombined.depthWindowCount() + region.depthWindowCount());
    }

    public void setTumorCopyNumber(@NotNull final CopyNumberMethod method, double copyNumber)
    {
        setCopyNumberMethod(method);
        mCombined.setTumorCopyNumber(copyNumber);
    }

    public void setInferredTumorBAF(double baf)
    {
        mInferredBAF = true;
        mCombined.setTumorBAF(baf);
        mCombined.setBafCount(0);
        mCombined.setObservedBAF(0);
    }

    @Override
    public void resetDepthWindowCount()
    {
        mCombined.setDepthWindowCount(0);
    }

    private double weightedAverage(long currentWeight, double currentValue, long newWeight, double newValue)
    {
        if(Doubles.isZero(currentValue))
        {
            return newValue;
        }

        long totalWeight = currentWeight + newWeight;
        return (currentWeight * currentValue + newWeight * newValue) / totalWeight;
    }

    private void clearBAFValues()
    {
        mCombined.setBafCount(0);
    }

}
