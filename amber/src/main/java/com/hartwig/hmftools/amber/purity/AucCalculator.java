package com.hartwig.hmftools.amber.purity;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

import com.hartwig.hmftools.amber.PositionEvidence;
import com.hartwig.hmftools.amber.contamination.CategoryEvidence;
import com.hartwig.hmftools.amber.contamination.CategoryEvidenceIntegral;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.segmentation.ChrArm;
import com.hartwig.hmftools.common.segmentation.ChrArmLocator;

public class AucCalculator<S, T extends Comparable<T>>
{
    public static AucCalculator<PositionEvidence, ChrArm> chrArmAucCalculator(RefGenomeVersion genomeVersion,
            List<PositionEvidence> hetPoints)
    {
        ChrArmLocator locator = ChrArmLocator.defaultLocator(genomeVersion);
        Function<PositionEvidence, ChrArm> classifier = evidence -> locator.map(evidence.chromosome(), evidence.position());
        return new AucCalculator<>(classifier, hetPoints);
    }

    public static AucCalculator<PositionEvidence, CanonicalSnvType> mutationTypeAucCalculator(List<PositionEvidence> hetPoints)
    {
        Function<PositionEvidence, CanonicalSnvType> classifier = evidence ->
        {
            if(evidence.vaf() < 0.5)
            {
                return CanonicalSnvType.type(evidence.Ref, evidence.Alt);
            }
            else
            {
                return CanonicalSnvType.type(evidence.Alt, evidence.Ref);
            }
        };
        return new AucCalculator<>(classifier, hetPoints);
    }

    private final Map<T, Integer> BaselineValues;
    private final Function<S, T> Classifier;

    public AucCalculator(Function<S, T> classifier, Collection<S> baselines)
    {
        Classifier = classifier;
        BaselineValues = getCounts(baselines);
    }

    public double calculateAuc(Collection<S> samples)
    {
        Map<T, Integer> counts = getCounts(samples);
        Set<T> categories = new HashSet<>(counts.keySet());
        categories.addAll(BaselineValues.keySet());

        Set<CategoryEvidence<T>> evidenceValues = new HashSet<>();
        for(T category : categories)
        {
            CategoryEvidence<T> evidence = new CategoryEvidence<>(category);
            int baselineValue = BaselineValues.getOrDefault(category, 0);
            int sampleValue = counts.getOrDefault(category, 0);
            evidence.set(baselineValue, sampleValue);
            evidenceValues.add(evidence);
        }
        CategoryEvidenceIntegral<T> integral = new CategoryEvidenceIntegral<>(evidenceValues);

        final int totalHits = integral.totalHits();
        final int totalPoints = integral.totalPoints();
        double possibleMax = 0.5 * totalHits * totalPoints;
        double cost = integral.value() / possibleMax;
        return cost;
    }

    private Map<T, Integer> getCounts(Collection<S> samples)
    {
        Map<T, Integer> counts = new HashMap<>();
        for(S s : samples)
        {
            final T category = Classifier.apply(s);
            int baselineValue = counts.getOrDefault(category, 0);
            counts.put(category, baselineValue + 1);
        }
        return counts;
    }
}
