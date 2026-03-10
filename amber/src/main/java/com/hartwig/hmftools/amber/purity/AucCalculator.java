package com.hartwig.hmftools.amber.purity;

import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.Function;

public class AucCalculator<S, T extends Comparable<T>>
{
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
        return integral.value() / possibleMax;
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
