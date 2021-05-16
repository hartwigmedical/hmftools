package com.hartwig.hmftools.linx.visualiser.data;

import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.utils.Doubles;
import com.hartwig.hmftools.linx.visualiser.file.VisCopyNumberFile;

import org.jetbrains.annotations.NotNull;

public class VisCopyNumbers
{
    private static final int MAX_EXTRA_DISTANCE = 1000;
    private static final double MIN_EXTRA_DISTANCE_PERCENT = 0.1;

    // From DB: select sampleId, chromosome, start, end, copyNumber, baf from copyNumber where sampleId = 'xxxxx';

    @NotNull
    public static List<CopyNumberAlteration> copyNumbers(@NotNull final List<CopyNumberAlteration> alterations,
            @NotNull final List<GenomeRegion> span)
    {
        final List<CopyNumberAlteration> result = Lists.newArrayList();

        for (int i = 0; i < alterations.size(); i++)
        {
            CopyNumberAlteration alteration = alterations.get(i);
            final String contig = alteration.chromosome();
            final List<GenomeRegion> chromosomeSegments =
                    span.stream().filter(x -> x.chromosome().equals(contig)).collect(Collectors.toList());
            if (!chromosomeSegments.isEmpty())
            {
                long minTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::start).min().orElse(0);
                long maxTrackPosition = chromosomeSegments.stream().mapToLong(GenomeRegion::end).max().orElse(0);
                long chromosomeDistance = maxTrackPosition - minTrackPosition;
                long additional = Math.max(1, Math.min(MAX_EXTRA_DISTANCE, Math.round(MIN_EXTRA_DISTANCE_PERCENT * chromosomeDistance)));
                minTrackPosition = minTrackPosition - additional;
                maxTrackPosition = maxTrackPosition + additional;

                if (alteration.end() >= minTrackPosition && alteration.start() <= maxTrackPosition)
                {
                    boolean isStartDecreasing = i > 0 && lessThan(alteration, alterations.get(i - 1));
                    long startPosition = isStartDecreasing ? alteration.start() - 1 : alteration.start();

                    boolean isEndIncreasing = i < alterations.size() - 1 && lessThan(alteration, alterations.get(i + 1));
                    long endPosition = isEndIncreasing ? alteration.end() + 1 : alteration.end();

                    result.add(ImmutableCopyNumberAlteration.builder()
                            .from(alteration)
                            .start(Math.max(minTrackPosition, startPosition))
                            .end(Math.min(maxTrackPosition, endPosition))
                            .truncated(minTrackPosition > startPosition || maxTrackPosition < endPosition)
                            .build());
                }
            }
        }

        return result;
    }

    private static boolean lessThan(@NotNull final CopyNumberAlteration first, @NotNull final CopyNumberAlteration second)
    {
        return first.chromosome().equals(second.chromosome()) && Doubles.lessThan(first.copyNumber(), second.copyNumber());
    }

    @NotNull
    public static List<CopyNumberAlteration> read(@NotNull final String fileName) throws IOException
    {
        return VisCopyNumberFile.read(fileName).stream().map(VisCopyNumbers::fromVis).collect(Collectors.toList());
    }

    @NotNull
    private static CopyNumberAlteration fromVis(@NotNull final VisCopyNumberFile visCopyNumberFile)
    {
        return ImmutableCopyNumberAlteration.builder()
                .sampleId(visCopyNumberFile.SampleId)
                .chromosome(visCopyNumberFile.Chromosome)
                .start(visCopyNumberFile.Start)
                .end(visCopyNumberFile.End)
                .copyNumber(visCopyNumberFile.CopyNumber)
                .baf(visCopyNumberFile.BAF)
                .truncated(false)
                .build();
    }

}
