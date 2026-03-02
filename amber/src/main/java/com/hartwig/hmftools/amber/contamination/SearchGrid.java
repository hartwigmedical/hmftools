package com.hartwig.hmftools.amber.contamination;

import static com.hartwig.hmftools.amber.AmberConfig.AMB_LOGGER;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Stream;

import com.google.common.base.Preconditions;

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

    public interface Calculator
    {
        ValueScore calculate(double value);
    }

    public List<Double> searchValues()
    {
        final double start = 0.005;
        final double end = 0.35;
        final double stepRatio = 1.05;
        double step = 0.001;
        List<Double> result = new ArrayList<>(60);
        double current = start;
        while(current <= end + 0.0001)
        {
            result.add(Math.round(current * 1000.0) / 1000.0);
            current += step;
            step *= stepRatio;
        }
        return result;
    }

    public ValueScore findBestValue(final Calculator calculator)
    {
        List<ValueScore> scores = new ArrayList<>();
        for(double value : searchValues())
        {
            ValueScore valueScore = calculator.calculate(value);
            scores.add(valueScore);
            //            System.out.println(value + "\t" + valueScore.score());
        }
        final List<ValueScore> sorted = scores.stream().sorted().toList();
        if(sorted.isEmpty())
        {
            return null;
        }
        AMB_LOGGER.debug("Best contamination score is {}, worst is {}", sorted.get(0), sorted.get(sorted.size() - 1));
        return sorted.get(0);
    }
}

