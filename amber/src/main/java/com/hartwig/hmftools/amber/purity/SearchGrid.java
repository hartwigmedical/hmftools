package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;

public class SearchGrid
{
    public record ValueScore(double value, double score) implements Comparable<ValueScore>
    {
        @Override
        public int compareTo(@NotNull final ValueScore o)
        {
            int byScore = Double.compare(score(), o.score());
            if(byScore == 0)
            {
                return Double.compare(value(), o.value());
            }
            return byScore;
        }
    }

    public List<Pair<Double, Double>> searchValuesAndSteps()
    {
        final double start = 0.005;
        final double end = 0.37;
        final double stepRatio = 1.05;
        double step = 0.001;
        List<Pair<Double, Double>> result = new ArrayList<>(60);
        double current = start;
        while(current <= end + 0.0001)
        {
            final double currentValue = Math.round(current * 1000.0) / 1000.0;
            current += step;
            step *= stepRatio;
            result.add(Pair.of(currentValue, step));
        }
        return result;
    }
}

