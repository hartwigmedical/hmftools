package com.hartwig.hmftools.common.gene;

import java.util.Collection;
import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.hmfslicer.HmfExonRegion;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;

import org.jetbrains.annotations.NotNull;

public class CanonicalTranscriptFactory {

    @NotNull
    public static List<CanonicalTranscript> create(@NotNull final Collection<HmfGenomeRegion> regions) {

        final List<CanonicalTranscript> transcripts = Lists.newArrayList();
        for (final HmfGenomeRegion region : regions) {

            final CanonicalTranscript transcript = create(region);

            transcripts.add(transcript);
        }

        return transcripts;
    }

    @VisibleForTesting
    @NotNull
    static CanonicalTranscript create(final HmfGenomeRegion region) {
        long codingStart = region.codingStart();
        long codingEnd = region.codingEnd();
        long exonStart = 0;
        long exonEnd = 0;
        long exonBases = 0;
        int codingBases = 0;
        int codingExons = 0;

        for (final HmfExonRegion exon : region.exome()) {
            if (exonStart == 0) {
                exonStart = exon.start();
            }

            if (codingStart <= exon.end() && codingEnd >= exon.start()) {
                codingExons++;
                codingBases += Math.min(exon.end(), codingEnd) - Math.max(exon.start(), codingStart) + 1;
            }

            exonBases += exon.bases();
            exonEnd = exon.end();
        }

        return ImmutableCanonicalTranscript.builder()
                .from(region)
                .geneID(region.geneID())
                .geneStart(region.geneStart())
                .geneEnd(region.geneEnd())
                .gene(region.gene())
                .exons(region.exome().size())
                .codingExons(codingExons)
                .exonStart(exonStart)
                .exonEnd(exonEnd)
                .exonBases(exonBases)
                .codingStart(codingStart)
                .codingEnd(codingEnd)
                .codingBases(Math.max(0, codingBases - 3))
                .build();
    }

}
