package com.hartwig.hmftools.common.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Modifiable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public abstract class HmfTranscriptRegion implements TranscriptRegion {

    @NotNull
    public abstract String geneID();

    public abstract long geneStart();

    public abstract long geneEnd();

    public abstract long codingStart();

    public abstract long codingEnd();

    @NotNull
    public abstract Strand strand();

    @NotNull
    public abstract List<Integer> entrezId();

    @NotNull
    public abstract List<HmfExonRegion> exome();

    @Value.Derived
    @Nullable
    public HmfExonRegion exonByIndex(int index) {
        int effectiveIndex = index - 1;
        List<HmfExonRegion> strandSortedExome = strandSortedExome();

        if (effectiveIndex >= 0 && effectiveIndex < strandSortedExome.size()) {
            return strandSortedExome.get(effectiveIndex);
        }

        return null;
    }

    @Value.Derived
    @Nullable
    public List<GenomeRegion> codonByIndex(int index) {
        return codonRangeByIndex(index, index);
    }

    @Value.Derived
    @Nullable
    public List<GenomeRegion> codonRangeByIndex(int startCodon, int endCodon) {
        if (startCodon < 1 || endCodon < 1) {
            // KODU: Enforce 1-based codons.
            return null;
        }

        if (codingStart() == 0 || codingEnd() == 0) {
            // KODU: Only coding transcripts have codons.
            return null;
        }

        List<GenomeRegion> codonRegions = Lists.newArrayList();
        int effectiveStartBase = 1 + (startCodon - 1) * 3;
        int effectiveEndBase = 3 + (endCodon - 1) * 3;

        int basesCovered = 0;
        Long startPosition = null;
        Long endPosition = null;
        for (HmfExonRegion exon : strandSortedExome()) {
            long exonCodingStart = Math.max(exon.start(), codingStart());
            long exonCodingEnd = Math.min(exon.end(), codingEnd());
            if (exonCodingStart > exonCodingEnd) {
                continue;
            }

            if (basesCovered + exon.bases() >= effectiveStartBase && startPosition == null) {
                startPosition = strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveStartBase - basesCovered - 1
                        : exonCodingEnd - effectiveStartBase + basesCovered + 1;
            }

            if (basesCovered + exon.bases() >= effectiveEndBase && endPosition == null) {
                endPosition = strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveEndBase - basesCovered - 1
                        : exonCodingEnd - effectiveEndBase + basesCovered + 1;
            }

            if (startPosition != null) {
                if (endPosition == null) {
                    codonRegions.add(ImmutableGenomeRegionImpl.builder()
                            .chromosome(chromosome())
                            .start(strand() == Strand.FORWARD ? startPosition : exonCodingStart)
                            .end(strand() == Strand.FORWARD ? exonCodingEnd : startPosition)
                            .build());
                } else if (codonRegions.size() == 1) {
                    codonRegions.add(ImmutableGenomeRegionImpl.builder()
                            .chromosome(chromosome())
                            .start(strand() == Strand.FORWARD ? exonCodingStart : endPosition)
                            .end(strand() == Strand.FORWARD ? endPosition : exonCodingEnd)
                            .build());
                    Collections.sort(codonRegions);
                    return codonRegions;
                } else {
                    codonRegions.add(ImmutableGenomeRegionImpl.builder()
                            .chromosome(chromosome())
                            .start(strand() == Strand.FORWARD ? startPosition :endPosition)
                            .end(strand() == Strand.FORWARD ? endPosition : startPosition)
                            .build());
                    return codonRegions;
                }
            }

            basesCovered += (1 + exonCodingEnd - exonCodingStart);
        }

        return null;
    }

    @Value.Derived
    @NotNull
    public List<HmfExonRegion> strandSortedExome() {
        return strand() == Strand.FORWARD ? exome() : Lists.reverse(exome());
    }
}
