package com.hartwig.hmftools.panelbuilder;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.panelbuilder.RegionUtils.regionOverlapsOrAdjacent;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.region.BaseRegion;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;

// Shared gene/transcript helpers used by the gene probe generators.
public class GeneUtils
{
    // A contiguous genome region formed by merging overlapping or adjacent exons across one or more transcripts, with the merged coding span
    // (if any part of the region is coding).
    public static class MergedExonRegion
    {
        public BaseRegion Region;
        public boolean IsCoding = false;
        public int CodingStart = Integer.MAX_VALUE;
        public int CodingEnd = Integer.MIN_VALUE;
        public List<ExonData> Exons = new ArrayList<>();
    }

    // Merges the exons of the given transcripts into a non-overlapping, ascending list of regions (overlapping or adjacent exons combined).
    public static List<MergedExonRegion> mergeExons(final List<TranscriptData> transcripts)
    {
        List<MergedExonRegion> mergedExons = new ArrayList<>();
        // Standard region merging algorithm...
        List<ImmutablePair<TranscriptData, ExonData>> sortedExons = transcripts.stream()
                .flatMap(trans -> trans.exons().stream().map(exon -> new ImmutablePair<>(trans, exon)))
                .sorted(Comparator.comparing(pair -> pair.getRight().Start))
                .toList();
        for(Pair<TranscriptData, ExonData> pair : sortedExons)
        {
            TranscriptData transcript = pair.getLeft();
            ExonData exon = pair.getRight();

            BaseRegion codingRegion = transcript.nonCoding() ? null : new BaseRegion(transcript.CodingStart, transcript.CodingEnd);
            BaseRegion exonRegion = new BaseRegion(exon.Start, exon.End);
            boolean isCoding = codingRegion != null && codingRegion.overlaps(exonRegion);

            mergedExons.stream()
                    .filter(region -> regionOverlapsOrAdjacent(region.Region, exonRegion))
                    .findFirst()
                    .ifPresentOrElse(
                            merged ->
                            {
                                if(exon.End > merged.Region.end())
                                {
                                    // Don't mutate in place because we borrowed the object from the exon data.
                                    merged.Region = new BaseRegion(merged.Region.start(), exon.End);
                                }
                                if(isCoding)
                                {
                                    merged.IsCoding = true;
                                    merged.CodingStart = min(merged.CodingStart, codingRegion.start());
                                    merged.CodingEnd = max(merged.CodingEnd, codingRegion.end());
                                }
                                merged.Exons.add(exon);
                            },
                            () ->
                            {
                                MergedExonRegion merged = new MergedExonRegion();
                                merged.Region = exonRegion;
                                merged.IsCoding = isCoding;
                                if(isCoding)
                                {
                                    merged.CodingStart = codingRegion.start();
                                    merged.CodingEnd = codingRegion.end();
                                }
                                merged.Exons.add(exon);
                                mergedExons.add(merged);
                            });
        }
        return mergedExons;
    }
}
