package com.hartwig.hmftools.amber.purity;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.amber.AmberConstants;

import org.apache.commons.lang3.tuple.Pair;

public class SearchGrid
{

    public record ValueScore(double value, double score) implements Comparable<ValueScore>
    {
        @Override
        public int compareTo(final ValueScore o)
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
        double step = AmberConstants.PEAK_SEARCH_INITIAL_STEP;
        List<Pair<Double, Double>> result = new ArrayList<>(60);
        double current = AmberConstants.PEAK_SEARCH_START;
        while(current <= AmberConstants.PEAK_SEARCH_END + AmberConstants.PEAK_SEARCH_OVERSHOOT)
        {
            final double currentValue = Math.round(current * 1000.0) / 1000.0;
            current += step;
            step *= AmberConstants.PEAK_SEARCH_STEP_RATIO;
            result.add(Pair.of(currentValue, step));
        }
        return result;
    }
}

