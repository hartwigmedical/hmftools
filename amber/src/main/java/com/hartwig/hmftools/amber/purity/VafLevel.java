package com.hartwig.hmftools.amber.purity;

import static java.util.stream.Collectors.toSet;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.VafReading;
import com.hartwig.hmftools.amber.contamination.CategoryEvidence;
import com.hartwig.hmftools.amber.contamination.PerArmVafConsistencyChecker;
import com.hartwig.hmftools.amber.contamination.VafPredicate;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

import org.apache.commons.math3.distribution.BinomialDistribution;

public class VafLevel
{
    private final double Level;
    private final List<PositionEvidence> PointsTested = new ArrayList<>();
    private final List<PositionEvidence> PointsInHomozygousBand = new ArrayList<>();
    private final List<PositionEvidence> PointsInHeterozygousBand = new ArrayList<>();

    public VafLevel(double Level)
    {
        this.Level = Level;
    }

    public boolean hasSufficientDepthForEventDetection(PositionEvidence evidence)
    {
        BinomialDistribution binomial = new BinomialDistribution(evidence.ReadDepth, Level);
        return binomial.cumulativeProbability(2) < 0.16;
    }

    public void test(PositionEvidence evidence)
    {
        PointsTested.add(evidence);

        int n = evidence.ReadDepth;
        int k = evidence.AltSupport;
        if(k > n / 2)
        {
            k = n - k;
        }
        BinomialDistribution binomialHomozygous = new BinomialDistribution(n, Level);
        double cdfHomozygous = binomialHomozygous.cumulativeProbability(k);
        if(cdfHomozygous > 0.25 && cdfHomozygous < 0.75)
        {
            PointsInHomozygousBand.add(evidence);
        }
        else
        {
            BinomialDistribution binomialHeterozygous = new BinomialDistribution(n, Level / 2);
            double cdfHeterozygous = binomialHeterozygous.cumulativeProbability(k);
            if(cdfHeterozygous > 0.25 && cdfHeterozygous < 0.75)
            {
                PointsInHeterozygousBand.add(evidence);
            }
        }
    }

    public double perArmConsistencyFactor(final ChrArmLocator chrArmLocator)
    {
        final Set<VafReading> capturedEvents = allCapturedPoints().stream().map(PositionEvidence::convertToVafReading).collect(toSet());
        VafPredicate classifier = capturedEvents::contains;
        PerArmVafConsistencyChecker checker = new PerArmVafConsistencyChecker(classifier, chrArmLocator);
        for(PositionEvidence evidence : PointsTested)
        {
            checker.offer(evidence.convertToVafReading());
        }
        return checker.unevenDistributionCost().unevenDistributionCost();
    }

    public double perMutationTypeConsistencyFactor()
    {
        return 0;
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

    private Set<PositionEvidence> allCapturedPoints()
    {
        Set<PositionEvidence> result = new HashSet<>();
        result.addAll(PointsInHomozygousBand);
        result.addAll(PointsInHeterozygousBand);
        return result;
    }
}
