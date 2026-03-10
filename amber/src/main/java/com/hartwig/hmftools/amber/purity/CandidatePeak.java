package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.amber.PositionEvidence;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class CandidatePeak
{
    private static final double LOWER_CDF_BOUND_FOR_CAPTURE = 0.16;
    private static final double UPPER_CDF_BOUND_FOR_CAPTURE = 0.84;
    private final double Level;
    private final double StepToNextLevel;
    private final List<PositionEvidence> PointsTested = new ArrayList<>();
    private final List<PositionEvidence> PointsInHomozygousBand = new ArrayList<>();
    private final List<PositionEvidence> PointsInHeterozygousBand = new ArrayList<>();

    public CandidatePeak(double Level)
    {
        this.Level = Level;
        StepToNextLevel = 0.00001;
    }

    public CandidatePeak(double Level, double StepToNextLevel)
    {
        this.Level = Level;
        this.StepToNextLevel = StepToNextLevel;
    }

    public boolean hasSufficientDepthForEventDetection(PositionEvidence evidence)
    {
        // We want het points to be above the depth-2 noise floor.
        BinomialDistribution binomial = new BinomialDistribution(evidence.RefSupport + evidence.AltSupport, Level / 2);
        return binomial.cumulativeProbability(2) < 0.16;
    }

    public void test(PositionEvidence evidence)
    {
        PointsTested.add(evidence);

        int n = evidence.RefSupport + evidence.AltSupport;
        int k = evidence.AltSupport;
        if(k > n / 2)
        {
            k = n - k;
        }
        RangeStepPeak homPeak = new RangeStepPeak(n, k, Level, StepToNextLevel);
        final double hetLevel = Level / 2;
        RangeStepPeak hetPeak = new RangeStepPeak(n, k, hetLevel, StepToNextLevel / 2);
        final double symmetricVaf = evidence.symmetricVaf();
        if(homPeak.isCaptured(symmetricVaf))
        {
            if(hetPeak.isCaptured(symmetricVaf))
            {
                // Add it to the band to which it is closest.
                if(symmetricVaf > Level * 0.75)
                {
                    PointsInHomozygousBand.add(evidence);
                }
                else
                {
                    PointsInHeterozygousBand.add(evidence);
                }
            }
            else
            {
                PointsInHomozygousBand.add(evidence);
            }
        }
        else
        {
            if(hetPeak.isCaptured(symmetricVaf))
            {
                PointsInHeterozygousBand.add(evidence);
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
        return new HashSet<>(PointsInHomozygousBand);
    }

    public Set<PositionEvidence> heterozygousEvidencePoints()
    {
        return new HashSet<>(PointsInHeterozygousBand);
    }

    public int numberOfCapturedEvidencePoints()
    {
        return allCapturedPoints().size();
    }

    public double vaf()
    {
        return Level;
    }

    public double homozygousProportion()
    {
        return (double) PointsInHomozygousBand.size() / (allCapturedPoints().size());
    }

    @Override
    public String toString()
    {
        return String.format("VafLevel{vaf=%.2f, step=%.2f, tested: %d, homozygous: %d, heterozygous: %d}", Level, StepToNextLevel, PointsTested.size(), PointsInHomozygousBand.size(), PointsInHeterozygousBand.size());
    }

    public Set<PositionEvidence> allCapturedPoints()
    {
        Set<PositionEvidence> result = new HashSet<>();
        result.addAll(PointsInHomozygousBand);
        result.addAll(PointsInHeterozygousBand);
        return result;
    }
}
