package com.hartwig.hmftools.amber.blacklist;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeSet;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.genome.position.GenomePositionImpl;

public class AmberBlacklistStatistics
{
    private static final double LOWER_BOUND_OF_MIDDLE_BAND = 0.4;
    private static final double UPPER_BOUND_OF_MIDDLE_BAND = 0.6;
    private final Map<GenomePositionImpl, Integer> Counts = new HashMap<>();
    private final Map<GenomePositionImpl, Double> Totals = new HashMap<>();

    public void record(GenomePosition position, Double vaf)
    {
        final GenomePositionImpl gp = toImpl(position);
        Counts.put(gp, Counts.getOrDefault(gp, 0) + 1);
        Totals.put(gp, Totals.getOrDefault(gp, 0.0) + vaf);
    }

    public List<AmberBlacklistPoint> findSuspiciousPoints(int numberOfSamples)
    {
        SortedSet<AmberBlacklistPoint> result = new TreeSet<>();
        int minCount = Math.max(10, numberOfSamples / 5);
        for(Map.Entry<GenomePositionImpl, Integer> entry : Counts.entrySet())
        {
            final Integer pointCount = entry.getValue();
            if(pointCount > minCount)
            {
                double averageVaf = Totals.get(entry.getKey()) / pointCount;
                if(averageVaf < LOWER_BOUND_OF_MIDDLE_BAND || averageVaf > UPPER_BOUND_OF_MIDDLE_BAND)
                {
                    result.add(new AmberBlacklistPoint(entry.getKey().chromosome(), entry.getKey().position(), pointCount, averageVaf));
                }
            }
        }
        return result.stream().toList();
    }

    private static GenomePositionImpl toImpl(GenomePosition genomePosition)
    {
        return new GenomePositionImpl(genomePosition.chromosome(), genomePosition.position());
    }
}
