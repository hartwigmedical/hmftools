package com.hartwig.hmftools.linx.visualiser.data;

import static java.lang.Math.max;
import static java.lang.Math.min;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumber;

public class VisCopyNumbers
{
    private static final int MAX_EXTRA_DISTANCE = 1000;
    private static final double MIN_EXTRA_DISTANCE_PERCENT = 0.1;
    
    public static List<VisCopyNumber> copyNumbers(final List<VisCopyNumber> copyNumbers, final List<GenomeRegion> span)
    {
        final List<VisCopyNumber> results = Lists.newArrayList();

        for(int i = 0; i < copyNumbers.size(); i++)
        {
            VisCopyNumber alteration = copyNumbers.get(i);
            final String contig = alteration.chromosome();
            final List<GenomeRegion> chromosomeSegments =
                    span.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if(!chromosomeSegments.isEmpty())
            {
                int minTrackPosition = chromosomeSegments.stream().mapToInt(GenomeRegion::start).min().orElse(0);
                int maxTrackPosition = chromosomeSegments.stream().mapToInt(GenomeRegion::end).max().orElse(0);
                int chromosomeDistance = maxTrackPosition - minTrackPosition;
                int additional = (int) max(1, min(MAX_EXTRA_DISTANCE, Math.round(MIN_EXTRA_DISTANCE_PERCENT * chromosomeDistance)));
                minTrackPosition = minTrackPosition - additional;
                maxTrackPosition = maxTrackPosition + additional;

                if(alteration.end() >= minTrackPosition && alteration.start() <= maxTrackPosition)
                {
                    boolean isStartDecreasing = i > 0 && lessThan(alteration, copyNumbers.get(i - 1));
                    int startPosition = isStartDecreasing ? alteration.start() - 1 : alteration.start();

                    boolean isEndIncreasing = i < copyNumbers.size() - 1 && lessThan(alteration, copyNumbers.get(i + 1));
                    int endPosition = isEndIncreasing ? alteration.end() + 1 : alteration.end();

                    VisCopyNumber newCn = VisCopyNumber.from(alteration);
                    newCn.Start = max(minTrackPosition, startPosition);
                    newCn.End = min(maxTrackPosition, endPosition);
                    newCn.Truncated = minTrackPosition > startPosition || maxTrackPosition < endPosition;

                    results.add(newCn);
                }
            }
        }

        return results;
    }

    private static boolean lessThan(final VisCopyNumber first, final VisCopyNumber second)
    {
        return first.chromosome().equals(second.chromosome()) && Doubles.lessThan(first.CopyNumber, second.CopyNumber);
    }

    public static List<VisCopyNumber> read(final String fileName) throws IOException
    {
        return VisCopyNumber.read(fileName);
    }

}
