package com.hartwig.hmftools.common.genome.genepanel;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.genome.bed.NamedBed;
import com.hartwig.hmftools.common.genome.bed.NamedBedBuilder;
import com.hartwig.hmftools.common.genome.bed.NamedBedFactory;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;
import com.hartwig.hmftools.common.genome.region.GenomeRegionsBuilder;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class HmfExonPanelBed {

    private static final Logger LOGGER = LogManager.getLogger(HmfExonPanelBed.class);

    private HmfExonPanelBed() {
    }

    @NotNull
    public static List<GenomeRegion> createUnnamedCodingRegions(boolean includeUTR, @NotNull final Set<String> actionableGenes,
            @NotNull final List<HmfTranscriptRegion> regions) {
        GenomeRegionsBuilder builder = new GenomeRegionsBuilder();
        regions.stream()
                .distinct()
                .filter(transcript -> transcript.codingStart() != 0 && actionableGenes.contains(transcript.geneName()))
                .flatMap(transcript -> NamedBedFactory.codingRegions(includeUTR, transcript).stream())
                .sorted()
                .forEach(builder::addRegion);
        return builder.build();
    }

    @NotNull
    public static List<NamedBed> createNamedCodingRegions(boolean includeUTR, @NotNull final Set<String> actionableGenes,
            @NotNull final List<HmfTranscriptRegion> regions) {
        final List<NamedBed> namedExons = regions.stream()
                .distinct()
                .filter(transcript -> transcript.codingStart() != 0 && actionableGenes.contains(transcript.geneName()))
                .flatMap(transcript -> NamedBedFactory.codingRegions(includeUTR, transcript).stream())
                .sorted()
                .collect(Collectors.toList());

        NamedBedBuilder builder = new NamedBedBuilder();
        for (NamedBed exon : namedExons) {
            if (!builder.addBed(exon)) {
                LOGGER.warn("Unable to add exon {}", exon);
            }
        }

        List<NamedBed> result = builder.build();

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
            return true;
        }

        for (int i = 0; i < sortedInput.size(); i++) {
            GenomeRegion first = sortedInput.get(i);
            GenomeRegion second = reduced.get(i);
            if (!first.equals(second)) {
                return true;
            }
        }

        return false;
    }
}
