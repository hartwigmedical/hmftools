package com.hartwig.hmftools.common.genome.bed;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.jetbrains.annotations.NotNull;

public final class NamedBedFactory
{
    private static final int SPLICE_SIZE = 10;

    private NamedBedFactory()
    {
    }

    @NotNull
    public static List<NamedBed> codingRegions(boolean includeUTR, @NotNull final HmfTranscriptRegion transcript)
    {
        final int startPosition = includeUTR ? transcript.start() : transcript.codingStart();
        final int endPosition = includeUTR ? transcript.end() : transcript.codingEnd();

        final List<NamedBed> result = Lists.newArrayList();
        final GenomeRegionsBuilder regionBuilder = new GenomeRegionsBuilder();

        for(int i = 0; i < transcript.exons().size(); i++)
        {
            final HmfExonRegion exon = transcript.exons().get(i);
            int exonStart = i == 0 ? exon.start() : exon.start() - SPLICE_SIZE;
            int exonEnd = i == transcript.exons().size() - 1 ? exon.end() : exon.end() + SPLICE_SIZE;

            if(startPosition < exonEnd && endPosition > exonStart)
            {
                regionBuilder.addRegion(exon.chromosome(), Math.max(startPosition, exonStart), Math.min(endPosition, exonEnd));
            }
        }

        regionBuilder.build()
                .stream()
                .map(x -> ImmutableNamedBed.builder().from(x).name(transcript.geneName()).build())
                .forEach(result::add);

        Collections.sort(result);
        return result;
    }
}
