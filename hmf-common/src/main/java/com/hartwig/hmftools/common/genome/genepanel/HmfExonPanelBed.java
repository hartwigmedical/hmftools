package com.hartwig.hmftools.common.genome.genepanel;

import java.util.Collections;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.bed.ImmutableNamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfExonRegion;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.genome.region.Strand;

import org.jetbrains.annotations.NotNull;

public final class HmfExonPanelBed {

    @NotNull
    public static List<NamedBed> createRegions(@NotNull final Set<String> actionableGenes,
            @NotNull final List<HmfTranscriptRegion> regions) {
        final List<NamedBed> result = regions.stream()
                .distinct()
                .filter(transcript -> transcript.codingStart() != 0 && actionableGenes.contains(transcript.gene()))
                .flatMap(transcript -> codingExonsWithSpliceSites(transcript).stream())
                .sorted()
                .collect(Collectors.toList());

        if (containsOverlappingSegments(result)) {
            throw new IllegalArgumentException("Generating invalid bed!");
        }

        return result;
    }

    private static boolean containsOverlappingSegments(@NotNull List<NamedBed> regions) {
        final List<GenomeRegion> sortedInput =
                regions.stream().map(x -> GenomeRegions.create(x.chromosome(), x.start(), x.end())).sorted().collect(Collectors.toList());

        final GenomeRegionsBuilder builder = new GenomeRegionsBuilder(1);
        sortedInput.forEach(builder::addRegion);

        final List<GenomeRegion> reduced = builder.build();
        if (sortedInput.size() != reduced.size()) {
            return false;
        }

        for (int i = 0; i < sortedInput.size(); i++) {
            GenomeRegion first = sortedInput.get(0);
            GenomeRegion second = reduced.get(0);
            if (!first.equals(second)) {
                return true;
            }
        }

        return false;

    }

    @NotNull
    static List<NamedBed> codingExonsWithSpliceSites(@NotNull final HmfTranscriptRegion transcript) {

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

            if (transcript.codingStart() < exon.end() && transcript.codingEnd() > exon.start()) {
                regionBuilder.addRegion(exon.chromosome(),
                        Math.max(transcript.codingStart(), exon.start()),
                        Math.min(transcript.codingEnd(), exon.end()));
            }
        }

        regionBuilder.build().stream().map(x -> ImmutableNamedBed.builder().from(x).name(transcript.gene()).build()).forEach(result::add);

        Collections.sort(result);
        return result;
    }
}
