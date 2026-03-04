package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.contamination.PerClassVafConsistencyChecker;
import com.hartwig.hmftools.amber.contamination.VafClassifier;
import com.hartwig.hmftools.amber.contamination.VafPredicate;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VafLevel
{
    private static final double LOWER_CDF_BOUND_FOR_CAPTURE = 0.16;
    private static final double UPPER_CDF_BOUND_FOR_CAPTURE = 0.84;
    private final double Level;
    private final double StepToNextLevel;
    private final List<PositionEvidence> PointsTested = new ArrayList<>();
    private final List<PositionEvidence> PointsInHomozygousBand = new ArrayList<>();
    private final List<PositionEvidence> PointsInHeterozygousBand = new ArrayList<>();

    public VafLevel(double Level)
    {
        this.Level = Level;
        StepToNextLevel = 0.00001;
    }

    public VafLevel(double Level, double StepToNextLevel)
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
        BinomialDistribution binomialHomozygous = new BinomialDistribution(n, Level);
        double cdfHomozygous = binomialHomozygous.cumulativeProbability(k);
        final boolean capturedByCdfHom = cdfHomozygous > LOWER_CDF_BOUND_FOR_CAPTURE && cdfHomozygous < UPPER_CDF_BOUND_FOR_CAPTURE;
        final boolean capturedByStepHom = Math.abs(evidence.symmetricVaf() - Level) < StepToNextLevel;
        if(capturedByCdfHom || capturedByStepHom)
        {
            PointsInHomozygousBand.add(evidence);
        }
        else
        {
            BinomialDistribution binomialHeterozygous = new BinomialDistribution(n, Level / 2);
            double cdfHeterozygous = binomialHeterozygous.cumulativeProbability(k);
            final boolean capturedByCdfHet = cdfHeterozygous > LOWER_CDF_BOUND_FOR_CAPTURE && cdfHeterozygous < UPPER_CDF_BOUND_FOR_CAPTURE;
            final boolean capturedByStepHet = Math.abs(evidence.symmetricVaf() - Level / 2) < StepToNextLevel / 2;
            if(capturedByCdfHet || capturedByStepHet)
            {
                PointsInHeterozygousBand.add(evidence);
            }
        }
    }

    public double perArmConsistencyFactor(final ChrArmLocator chrArmLocator)
    {
        VafClassifier<PositionEvidence, ChrArm> chrArmClassifier = VafClassifier.chrArmPositionEvidenceClassifier(chrArmLocator);
        return calculateAucForClassifier(chrArmClassifier);
    }

    public double perMutationTypeConsistencyFactor()
    {
        VafClassifier<PositionEvidence, CanonicalSnvType> mutationTypeClassifier = VafClassifier.mutationTypeClassifier();
        return calculateAucForClassifier(mutationTypeClassifier);

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

    private <T extends Comparable<T>> double calculateAucForClassifier(final VafClassifier<PositionEvidence, T> chrArmClassifier)
    {
        VafPredicate<PositionEvidence> classifier = allCapturedPoints()::contains;
        PerClassVafConsistencyChecker<PositionEvidence, T> checker = new PerClassVafConsistencyChecker<>(classifier, chrArmClassifier);
        for(PositionEvidence evidence : PointsTested)
        {
            checker.offer(evidence);
        }
        return checker.unevenDistributionCost().unevenDistributionCost();
    }
}
