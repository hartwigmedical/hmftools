package com.hartwig.hmftools.common.genome.bed;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;

public class NamedBedFactory {

    @NotNull
    public static List<NamedBed> codingRegions(boolean includeUTR, @NotNull final HmfTranscriptRegion transcript) {
        final long startPosition = includeUTR ? transcript.start() : transcript.codingStart();
        final long endPosition = includeUTR ? transcript.end() : transcript.codingEnd();

        final List<NamedBed> result = Lists.newArrayList();
        final GenomeRegionsBuilder regionBuilder = new GenomeRegionsBuilder();

        boolean forward = transcript.strand() == Strand.FORWARD;
        for (int i = 0; i < transcript.exome().size(); i++) {

            // Splice sites (+1,+2,+5, -2,-1)
            final HmfExonRegion exon = transcript.exome().get(i);
            if (i != 0) {
                regionBuilder.addRegion(exon.chromosome(), exon.start() - 2, exon.start() - 1);
                if (!forward) {
                    regionBuilder.addPosition(exon.chromosome(), exon.start() - 5);
                }
            }
            if (i != transcript.exome().size() - 1) {
                regionBuilder.addRegion(exon.chromosome(), exon.end() + 1, exon.end() + 2);
                if (forward) {
                    regionBuilder.addPosition(exon.chromosome(), exon.end() + 5);
                }
            }

            if (startPosition < exon.end() && endPosition > exon.start()) {
                regionBuilder.addRegion(exon.chromosome(), Math.max(startPosition, exon.start()), Math.min(endPosition, exon.end()));
            }
        }

        regionBuilder.build().stream().map(x -> ImmutableNamedBed.builder().from(x).name(transcript.gene()).build()).forEach(result::add);

        Collections.sort(result);
        return result;
    }

}
