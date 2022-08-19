package com.hartwig.hmftools.common.genome.region;

import java.util.Collections;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;

import org.jetbrains.annotations.Nullable;

public class HmfTranscriptRegionUtils
{
    public static HmfTranscriptRegion fromTranscript(final GeneData geneData, final TranscriptData transData)
    {
        List<HmfExonRegion> exons = Lists.newArrayList();

        for(ExonData exon : transData.exons())
        {
            exons.add(ImmutableHmfExonRegion.builder()
                    .chromosome(geneData.Chromosome)
                    .exonRank(exon.Rank)
                    .start(exon.Start)
                    .end(exon.End)
                    .build());
        }

        return ImmutableHmfTranscriptRegion.builder()
                .geneId(geneData.GeneId)
                .geneName(geneData.GeneName)
                .chromosome(geneData.Chromosome)
                .strand(Strand.valueOf(geneData.Strand))
                .geneStart(geneData.GeneStart)
                .geneEnd(geneData.GeneEnd)
                .entrezId(Lists.newArrayList())
                .chromosomeBand(geneData.KaryotypeBand)
                .transName(transData.TransName)
                .isCanonical(transData.IsCanonical)
                .start(transData.TransStart)
                .end(transData.TransEnd)
                .codingStart(transData.CodingStart != null ? transData.CodingStart : -1)
                .codingEnd(transData.CodingEnd != null ? transData.CodingEnd : -1)
                .exons(exons)
                .build();
    }

    @Nullable
    public static List<GenomeRegion> codonByIndex(final HmfTranscriptRegion transcript, int index)
    {
        return codonRangeByRank(transcript, index, index);
    }

    @Nullable
    public static List<GenomeRegion> codonRangeByRank(final HmfTranscriptRegion transcript, int startCodon, int endCodon)
    {
        if(startCodon < 1 || endCodon < 1)
        {
            // Enforce 1-based codons.
            return null;
        }

        if(transcript.codingStart() == 0 || transcript.codingEnd() == 0)
        {
            // Only coding transcripts have codons.
            return null;
        }

        List<GenomeRegion> codonRegions = Lists.newArrayList();
        int effectiveStartBase = 1 + (startCodon - 1) * 3;
        int effectiveEndBase = 3 + (endCodon - 1) * 3;

        int basesCovered = 0;
        Integer startPosition = null;
        Integer endPosition = null;
        for(HmfExonRegion exon : transcript.strandSortedExome())
        {
            int exonCodingStart = Math.max(exon.start(), transcript.codingStart());
            int exonCodingEnd = Math.min(exon.end(), transcript.codingEnd());
            int exonBaseLength = exonCodingEnd - exonCodingStart + 1;

            if(exonBaseLength <= 0)
            {
                // Exon is entirely non-coding so can be skipped.
                continue;
            }

            if(basesCovered + exonBaseLength >= effectiveStartBase && startPosition == null)
            {
                startPosition = transcript.strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveStartBase - basesCovered - 1
                        : exonCodingEnd - effectiveStartBase + basesCovered + 1;
            }

            if(basesCovered + exonBaseLength >= effectiveEndBase && endPosition == null)
            {
                endPosition = transcript.strand() == Strand.FORWARD
                        ? exonCodingStart + effectiveEndBase - basesCovered - 1
                        : exonCodingEnd - effectiveEndBase + basesCovered + 1;
            }

            basesCovered += exonBaseLength;

            GenomeRegion region =
                    decideOnRangeToIncludeForExon(
                            transcript, startPosition, endPosition, exonCodingStart, exonCodingEnd,
                            codonRegions.size() > 0);

            if(region != null)
            {
                codonRegions.add(region);
            }

            if(startPosition != null && endPosition != null)
            {
                Collections.sort(codonRegions);
                return codonRegions;
            }
        }

        return null;
    }

    @Nullable
    private static GenomeRegion decideOnRangeToIncludeForExon(
            final HmfTranscriptRegion transcript, @Nullable Integer startPosition, @Nullable Integer endPosition, int exonCodingStart,
            int exonCodingEnd, boolean hasCodingRegionsDefinedAlready)
    {
        if(startPosition != null)
        {
            if(endPosition == null)
            {
                // Check to see if we need to include the entire exon we are considering.
                if(hasCodingRegionsDefinedAlready)
                {
                    return ImmutableGenomeRegionImpl.builder()
                            .chromosome(transcript.chromosome())
                            .start(exonCodingStart)
                            .end(exonCodingEnd)
                            .build();
                }
                else
                {
                    return ImmutableGenomeRegionImpl.builder()
                            .chromosome(transcript.chromosome())
                            .start(transcript.strand() == Strand.FORWARD ? startPosition : exonCodingStart)
                            .end(transcript.strand() == Strand.FORWARD ? exonCodingEnd : startPosition)
                            .build();
                }
            }
            else if(hasCodingRegionsDefinedAlready)
            {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(transcript.chromosome())
                        .start(transcript.strand() == Strand.FORWARD ? exonCodingStart : endPosition)
                        .end(transcript.strand() == Strand.FORWARD ? endPosition : exonCodingEnd)
                        .build();
            }
            else
            {
                return ImmutableGenomeRegionImpl.builder()
                        .chromosome(transcript.chromosome())
                        .start(transcript.strand() == Strand.FORWARD ? startPosition : endPosition)
                        .end(transcript.strand() == Strand.FORWARD ? endPosition : startPosition)
                        .build();
            }
        }
        return null;
    }
}
