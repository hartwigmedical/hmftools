package com.hartwig.hmftools.common.genome.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;

import org.jetbrains.annotations.NotNull;

public class CodingRegions {

    @NotNull
    public static List<NamedBed> codingRegions(@NotNull final HmfTranscriptRegion transcript) {
        if (transcript.codingStart() == 0) {
            return Collections.emptyList();
        }

        final long startPosition = transcript.codingStart();
        final long endPosition = transcript.codingEnd();

        final List<NamedBed> result = Lists.newArrayList();

        for (int i = 0; i < transcript.exome().size(); i++) {
            final HmfExonRegion exon = transcript.exome().get(i);
            if (startPosition <= exon.end() && endPosition >= exon.start()) {
                result.add(ImmutableNamedBed.builder()
                        .chromosome(exon.chromosome())
                        .start(Math.max(startPosition, exon.start()))
                        .end(Math.min(endPosition, exon.end()))
                        .name(transcript.gene())
                        .build());
            }
        }

        return result;
    }

}
