package com.hartwig.hmftools.common.genome.region;

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

    public boolean isDonorMinusOne(long position) {
        return isDonorPlusOne(-1, position);
    }

    public boolean isDonorPlusFive(long position) {
        return isDonorPlusOne(4, position);
    }

    public boolean isAcceptorPlusThree(long position) {
        return isAcceptorPlusOne(2, position);
    }

    public boolean isDonorPlusOne(int offset, long position) {
        if (codingStart() == 0) {
            return false;
        }

        if (position < codingStart() || position > codingEnd()) {
            return false;
        }

        if (strand() == Strand.FORWARD) {
            for (int i = 0; i < exome().size() - 1; i++) {
                long donorSite = exome().get(i).end() + 1;
                if (position == donorSite + offset) {
                    return true;
                }
            }
        } else {
            for (int i = 1; i < exome().size(); i++) {
                long donorSite = exome().get(i).start() - 1;
                if (position == donorSite - offset) {
                    return true;
                }
            }
        }

        return false;
    }

    public boolean isAcceptorPlusOne(int offset, long position) {
        if (codingStart() == 0) {
            return false;
        }

        if (position < codingStart() || position > codingEnd()) {
            return false;
        }

        if (strand() == Strand.FORWARD) {
            for (int i = 1; i < exome().size(); i++) {
                long acceptorSide = exome().get(i).start() - 1;
                if (position == acceptorSide - offset) {
                    return true;
                }
            }
        } else {
            for (int i = 0; i < exome().size() - 1; i++) {
                long donorSite = exome().get(i).end() + 1;
                if (position == donorSite + offset) {
                    return true;
                }
            }
        }

        return false;
    }

    @Nullable
    public HmfExonRegion exonByIndex(int index) {
        int effectiveIndex = index - 1;
        List<HmfExonRegion> strandSortedExome = strandSortedExome();

        if (effectiveIndex >= 0 && effectiveIndex < strandSortedExome.size()) {
            return strandSortedExome.get(effectiveIndex);
        }

        return null;
    }

    @Nullable
    public List<GenomeRegion> codingRangeByGenomicCoordinates(int genomicStartPosition, int genomicEndPosition) {
        if (genomicStartPosition < 1 || genomicEndPosition < 1) {
            // Only allow a range to return if start and end are 1-based.
            return null;
        }

        if (codingStart() == 0 || codingEnd() == 0) {
            // Only coding transcripts have coding ranges.
            return null;
        }

        long effectiveStartPosition = Math.max(codingStart(), genomicStartPosition);
        long effectiveEndPosition = Math.min(codingEnd(), genomicEndPosition);

        List<GenomeRegion> codingRegions = Lists.newArrayList();
        for (HmfExonRegion exon : exome()) {
            if (exon.start() <= effectiveEndPosition && exon.end() >= effectiveStartPosition) {
                codingRegions.add(ImmutableGenomeRegionImpl.builder()
                        .chromosome(chromosome())
                        .start(Math.max(exon.start(), effectiveStartPosition))
                        .end(Math.min(exon.end(), effectiveEndPosition))
                        .build());
            }
        }

        return !codingRegions.isEmpty() ? codingRegions : null;
    }

    @Nullable
    public List<GenomeRegion> codonByIndex(int index) {
        return codonRangeByIndex(index, index);
    }

    @NotNull
    public List<GenomeRegion> codonRangeAtGenomicPosition(long position) {
        final List<GenomeRegion> codonRegions = Lists.newArrayList();
        if (position < codingStart() || position > codingEnd()) {
            return codonRegions;
        }

        int basesCovered = 0;
        for (int i = 0; i < exome().size(); i++) {
            final HmfExonRegion exon = exome().get(i);
            long exonCodingStart = Math.max(exon.start(), codingStart());
            long exonCodingEnd = Math.min(exon.end(), codingEnd());
            long exonBaseLength = exonCodingEnd - exonCodingStart + 1;

            if (exonBaseLength <= 0) {
                // Exon is entirely non-coding so can be skipped.
                continue;
            }

            if (position >= exonCodingStart && position <= exonCodingEnd) {
                long lookBack = (basesCovered + position - exonCodingStart) % 3;
                long lookForward = 2 - lookBack;

                // Do we need previous exon?
                if (position - lookBack < exon.start() && i > 0) {
                    final HmfExonRegion previous = exome().get(i - 1);
                    final long previousExonLookBack = lookBack + exon.start() - position - 1;
                    codonRegions.add(GenomeRegions.create(chromosome(),
                            Math.max(previous.start(), previous.end() - previousExonLookBack),
                            previous.end()));
                }

                // Current exon
                codonRegions.add(GenomeRegions.create(chromosome(),
                        Math.max(exon.start(), position - lookBack),
                        Math.min(exon.end(), position + lookForward)));

                // Do we need next exon?
                if (position + lookForward > exon.end() && i < exome().size() - 1) {
                    final HmfExonRegion next = exome().get(i + 1);
                    final long nextExonLookForward = lookForward - exon.end() + position - 1;
                    codonRegions.add(GenomeRegions.create(chromosome(),
                            next.start(),
                            Math.min(next.end(), next.start() + nextExonLookForward)));
                }

                return codonRegions;
            }

            if (exonCodingStart > position) {
                return codonRegions;
            }

            basesCovered += exonBaseLength;
        }

        return codonRegions;
    }

    @Nullable
    public List<GenomeRegion> codonRangeByIndex(int startCodon, int endCodon) {
        if (startCodon < 1 || endCodon < 1) {
            // Enforce 1-based codons.
            return null;
        }

        if (codingStart() == 0 || codingEnd() == 0) {
            // Only coding transcripts have codons.
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
            long exonBaseLength = exonCodingEnd - exonCodingStart + 1;

            if (exonBaseLength <= 0) {
                // Exon is entirely non-coding so can be skipped.
                continue;
            }

            if (basesCovered + exonBaseLength >= effectiveStartBase && startPosition == null) {
                startPosition = strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveStartBase - basesCovered - 1
                        : exonCodingEnd - effectiveStartBase + basesCovered + 1;
            }

            if (basesCovered + exonBaseLength >= effectiveEndBase && endPosition == null) {
                endPosition = strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveEndBase - basesCovered - 1
                        : exonCodingEnd - effectiveEndBase + basesCovered + 1;
            }

            basesCovered += exonBaseLength;

            GenomeRegion region =
                    decideOnRangeToIncludeForExon(startPosition, endPosition, exonCodingStart, exonCodingEnd, codonRegions.size() > 0);

            if (region != null) {
                codonRegions.add(region);
            }

            if (startPosition != null && endPosition != null) {
                Collections.sort(codonRegions);
                return codonRegions;
            }
        }

        return null;
    }

    @Nullable
    private GenomeRegion decideOnRangeToIncludeForExon(@Nullable Long startPosition, @Nullable Long endPosition, long exonCodingStart,
            long exonCodingEnd, boolean hasCodingRegionsDefinedAlready) {
        if (startPosition != null) {
            if (endPosition == null) {
                // Check to see if we need to include the entire exon we are considering.
                if (hasCodingRegionsDefinedAlready) {
                    return ImmutableGenomeRegionImpl.builder().chromosome(chromosome()).start(exonCodingStart).end(exonCodingEnd).build();
                } else {
                    return ImmutableGenomeRegionImpl.builder()
                            .chromosome(chromosome())
                            .start(strand() == Strand.FORWARD ? startPosition : exonCodingStart)
                            .end(strand() == Strand.FORWARD ? exonCodingEnd : startPosition)
                            .build();
                }
            } else if (hasCodingRegionsDefinedAlready) {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(chromosome())
                        .start(strand() == Strand.FORWARD ? exonCodingStart : endPosition)
                        .end(strand() == Strand.FORWARD ? endPosition : exonCodingEnd)
                        .build();
            } else {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(chromosome())
                        .start(strand() == Strand.FORWARD ? startPosition : endPosition)
                        .end(strand() == Strand.FORWARD ? endPosition : startPosition)
                        .build();
            }
        }
        return null;
    }

    @NotNull
    private List<HmfExonRegion> strandSortedExome() {
        return strand() == Strand.FORWARD ? exome() : Lists.reverse(exome());
    }
}
