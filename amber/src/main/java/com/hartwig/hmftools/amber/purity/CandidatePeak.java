package com.hartwig.hmftools.amber.purity;

import static com.hartwig.hmftools.amber.AmberConstants.LOWER_CDF_BOUND_FOR_CAPTURE;
import static com.hartwig.hmftools.amber.AmberConstants.UPPER_CDF_BOUND_FOR_CAPTURE;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class CandidatePeak
{
    private final double mLevel;
    private final double mStepToNextLevel;
    private final List<PositionEvidence> mPointsTested = new ArrayList<>();
    private final List<PositionEvidence> mPointsInHomozygousBand = new ArrayList<>();
    private final List<PositionEvidence> mPointsInHeterozygousBand = new ArrayList<>();

    public CandidatePeak(final double level)
    {
        this.mLevel = level;
        mStepToNextLevel = 0.00001;
    }

    public CandidatePeak(final double level, final double stepToNextLevel)
    {
        this.mLevel = level;
        this.mStepToNextLevel = stepToNextLevel;
    }

    public boolean hasSufficientDepthForEventDetection(final PositionEvidence evidence)
    {
        // We want het points to be above the depth-2 noise floor.
        BinomialDistribution binomial = new BinomialDistribution(evidence.RefSupport + evidence.AltSupport, mLevel / 2);
        return binomial.cumulativeProbability(2) < LOWER_CDF_BOUND_FOR_CAPTURE;
    }

    public void test(final PositionEvidence evidence)
    {
        mPointsTested.add(evidence);

        int n = evidence.RefSupport + evidence.AltSupport;
        int k = evidence.AltSupport;
        if(k > n / 2)
        {
            k = n - k;
        }
        RangeStepPeak homPeak = new RangeStepPeak(n, k, mLevel, mStepToNextLevel);
        final double hetLevel = mLevel / 2;
        RangeStepPeak hetPeak = new RangeStepPeak(n, k, hetLevel, mStepToNextLevel / 2);
        final double symmetricVaf = evidence.symmetricVaf();
        if(homPeak.isCaptured(symmetricVaf))
        {
            if(hetPeak.isCaptured(symmetricVaf))
            {
                // Add it to the band to which it is closest.
                if(symmetricVaf > mLevel * 0.75)
                {
                    mPointsInHomozygousBand.add(evidence);
                }
                else
                {
                    mPointsInHeterozygousBand.add(evidence);
                }
            }
            else
            {
                mPointsInHomozygousBand.add(evidence);
            }
        }
        else
        {
            if(hetPeak.isCaptured(symmetricVaf))
            {
                mPointsInHeterozygousBand.add(evidence);
            }
        }
    }

    record RangeStepPeak(int depth, int count, double vaf, double step)
    {
        public boolean isCaptured(double symmetricVaf)
        {
            BinomialDistribution binomialHomozygous = new BinomialDistribution(depth, vaf);
            double cdfHomozygous = binomialHomozygous.cumulativeProbability(count);
            final boolean capturedByCdfHom = cdfHomozygous > LOWER_CDF_BOUND_FOR_CAPTURE && cdfHomozygous < UPPER_CDF_BOUND_FOR_CAPTURE;
            final boolean capturedByStepHom = Math.abs(symmetricVaf - vaf) < step;
            return capturedByCdfHom || capturedByStepHom;
        }
    }

    public Set<PositionEvidence> homozygousEvidencePoints()
    {
        return new HashSet<>(mPointsInHomozygousBand);
    }

    public Set<PositionEvidence> heterozygousEvidencePoints()
    {
        return new HashSet<>(mPointsInHeterozygousBand);
    }

    public int numberOfCapturedEvidencePoints()
    {
        return allCapturedPoints().size();
    }

    public double vaf()
    {
        return mLevel;
    }

    public double homozygousProportion()
    {
        return (double) mPointsInHomozygousBand.size() / (allCapturedPoints().size());
    }

    @Override
    public String toString()
    {
        return String.format("VafLevel{vaf=%.2f, step=%.2f, tested: %d, homozygous: %d, heterozygous: %d}",
                mLevel, mStepToNextLevel, mPointsTested.size(), mPointsInHomozygousBand.size(), mPointsInHeterozygousBand.size());
    }

    public Set<PositionEvidence> allCapturedPoints()
    {
        Set<PositionEvidence> result = new HashSet<>();
        result.addAll(mPointsInHomozygousBand);
        result.addAll(mPointsInHeterozygousBand);
        return result;
    }
}
