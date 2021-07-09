package com.hartwig.hmftools.common.genome.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class HmfTranscriptRegionUtils
{
    @Nullable
    public static List<GenomeRegion> codonByIndex(final HmfTranscriptRegion transcript, int index) {
        return codonRangeByIndex(transcript, index, index);
    }

    @Nullable
    public static List<GenomeRegion> codonRangeByIndex(final HmfTranscriptRegion transcript, int startCodon, int endCodon) {
        if (startCodon < 1 || endCodon < 1) {
            // Enforce 1-based codons.
            return null;
        }

        if (transcript.codingStart() == 0 || transcript.codingEnd() == 0) {
            // Only coding transcripts have codons.
            return null;
        }

        List<GenomeRegion> codonRegions = Lists.newArrayList();
        int effectiveStartBase = 1 + (startCodon - 1) * 3;
        int effectiveEndBase = 3 + (endCodon - 1) * 3;

        int basesCovered = 0;
        Long startPosition = null;
        Long endPosition = null;
        for (HmfExonRegion exon : transcript.strandSortedExome()) {
            long exonCodingStart = Math.max(exon.start(), transcript.codingStart());
            long exonCodingEnd = Math.min(exon.end(), transcript.codingEnd());
            long exonBaseLength = exonCodingEnd - exonCodingStart + 1;

            if (exonBaseLength <= 0) {
                // Exon is entirely non-coding so can be skipped.
                continue;
            }

            if (basesCovered + exonBaseLength >= effectiveStartBase && startPosition == null) {
                startPosition = transcript.strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveStartBase - basesCovered - 1
                        : exonCodingEnd - effectiveStartBase + basesCovered + 1;
            }

            if (basesCovered + exonBaseLength >= effectiveEndBase && endPosition == null) {
                endPosition = transcript.strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveEndBase - basesCovered - 1
                        : exonCodingEnd - effectiveEndBase + basesCovered + 1;
            }

            basesCovered += exonBaseLength;

            GenomeRegion region =
                    decideOnRangeToIncludeForExon(
                            transcript, startPosition, endPosition, exonCodingStart, exonCodingEnd,
                            codonRegions.size() > 0);

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
    private static GenomeRegion decideOnRangeToIncludeForExon(
            final HmfTranscriptRegion transcript, @Nullable Long startPosition, @Nullable Long endPosition, long exonCodingStart,
            long exonCodingEnd, boolean hasCodingRegionsDefinedAlready) {
        if (startPosition != null) {
            if (endPosition == null) {
                // Check to see if we need to include the entire exon we are considering.
                if (hasCodingRegionsDefinedAlready) {
                    return ImmutableGenomeRegionImpl.builder().chromosome(transcript.chromosome()).start(exonCodingStart).end(exonCodingEnd).build();
                } else {
                    return ImmutableGenomeRegionImpl.builder()
                            .chromosome(transcript.chromosome())
                            .start(transcript.strand() == Strand.FORWARD ? startPosition : exonCodingStart)
                            .end(transcript.strand() == Strand.FORWARD ? exonCodingEnd : startPosition)
                            .build();
                }
            } else if (hasCodingRegionsDefinedAlready) {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(transcript.chromosome())
                        .start(transcript.strand() == Strand.FORWARD ? exonCodingStart : endPosition)
                        .end(transcript.strand() == Strand.FORWARD ? endPosition : exonCodingEnd)
                        .build();
            } else {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(transcript.chromosome())
                        .start(transcript.strand() == Strand.FORWARD ? startPosition : endPosition)
                        .end(transcript.strand() == Strand.FORWARD ? endPosition : startPosition)
                        .build();
            }
        }
        return null;
    }

    @NotNull
    public static List<GenomeRegion> codonRangeAtGenomicPosition(final HmfTranscriptRegion transcript, long position)
    {
        final List<GenomeRegion> codonRegions = Lists.newArrayList();

        if(position < transcript.codingStart() || position > transcript.codingEnd())
            return codonRegions;

        int basesCovered = 0;
        for(int i = 0; i < transcript.exome().size(); i++)
        {
            final HmfExonRegion exon = transcript.exome().get(i);
            long exonCodingStart = Math.max(exon.start(), transcript.codingStart());
            long exonCodingEnd = Math.min(exon.end(), transcript.codingEnd());
            long exonBaseLength = exonCodingEnd - exonCodingStart + 1;

            if(exonBaseLength <= 0)
            {
                // Exon is entirely non-coding so can be skipped.
                continue;
            }

            if(position >= exonCodingStart && position <= exonCodingEnd)
            {
                long lookBack = (basesCovered + position - exonCodingStart) % 3;
                long lookForward = 2 - lookBack;

                // Do we need previous exon?
                if(position - lookBack < exon.start() && i > 0)
                {
                    final HmfExonRegion previous = transcript.exome().get(i - 1);
                    final long previousExonLookBack = lookBack + exon.start() - position - 1;
                    codonRegions.add(GenomeRegions.create(transcript.chromosome(),
                            Math.max(previous.start(), previous.end() - previousExonLookBack),
                            previous.end()));
                }

                // Current exon
                codonRegions.add(GenomeRegions.create(transcript.chromosome(),
                        Math.max(exon.start(), position - lookBack),
                        Math.min(exon.end(), position + lookForward)));

                // Do we need next exon?
                if(position + lookForward > exon.end() && i < transcript.exome().size() - 1)
                {
                    final HmfExonRegion next = transcript.exome().get(i + 1);
                    final long nextExonLookForward = lookForward - exon.end() + position - 1;
                    codonRegions.add(GenomeRegions.create(transcript.chromosome(),
                            next.start(),
                            Math.min(next.end(), next.start() + nextExonLookForward)));
                }

                return codonRegions;
            }

            if(exonCodingStart > position)
            {
                return codonRegions;
            }

            basesCovered += exonBaseLength;
        }

        return codonRegions;
    }

    public static boolean isDonorMinusOne(final HmfTranscriptRegion transcript, long position) {
        return isDonorPlusOne(transcript, -1, position);
    }

    public static boolean isDonorPlusFive(final HmfTranscriptRegion transcript, long position) {
        return isDonorPlusOne(transcript, 4, position);
    }

    public static boolean isAcceptorPlusThree(final HmfTranscriptRegion transcript, long position) {
        return isAcceptorPlusOne(transcript, 2, position);
    }

    private static boolean isDonorPlusOne(final HmfTranscriptRegion transcript, int offset, long position) {
        if (transcript.codingStart() == 0) {
            return false;
        }

        if (position < transcript.codingStart() || position > transcript.codingEnd()) {
            return false;
        }

        if (transcript.strand() == Strand.FORWARD) {
            for (int i = 0; i < transcript.exome().size() - 1; i++) {
                long donorSite = transcript.exome().get(i).end() + 1;
                if (position == donorSite + offset) {
                    return true;
                }
            }
        } else {
            for (int i = 1; i < transcript.exome().size(); i++) {
                long donorSite = transcript.exome().get(i).start() - 1;
                if (position == donorSite - offset) {
                    return true;
                }
            }
        }

        return false;
    }

    private static boolean isAcceptorPlusOne(final HmfTranscriptRegion transcript, int offset, long position) {
        if (transcript.codingStart() == 0) {
            return false;
        }

        if (position < transcript.codingStart() || position > transcript.codingEnd()) {
            return false;
        }

        if (transcript.strand() == Strand.FORWARD) {
            for (int i = 1; i < transcript.exome().size(); i++) {
                long acceptorSide = transcript.exome().get(i).start() - 1;
                if (position == acceptorSide - offset) {
                    return true;
                }
            }
        } else {
            for (int i = 0; i < transcript.exome().size() - 1; i++) {
                long donorSite = transcript.exome().get(i).end() + 1;
                if (position == donorSite + offset) {
                    return true;
                }
            }
        }

        return false;
    }

}
