package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

// Per-chromosome merged exon spans, keyed in versioned-chromosome form (V38 -> chr1) to match lifted genomic coords.
public final class ExonRegionIndex
{
    private final Map<String, List<BaseRegion>> mExonsByChromosome;

    private ExonRegionIndex(final Map<String, List<BaseRegion>> exonsByChromosome)
    {
        mExonsByChromosome = exonsByChromosome;
    }

    public boolean contains(final String chromosome, final int pos)
    {
        List<BaseRegion> exons = mExonsByChromosome.get(chromosome);
        if(exons == null)
            return false;

        // merged spans are non-overlapping, so the last exon starting at or before pos is the only candidate.
        BaseRegion exon = exons.get(BaseRegion.binarySearch(pos, exons));
        return exon.containsPosition(pos);
    }

    // Built from the contig sidecar entries' exon spans (multi-exon transcripts only; single-exon transcripts
    // have no tx-contig, so are absent).
    public static ExonRegionIndex fromContigEntries(final List<ContigEntry> entries)
    {
        Map<String, List<BaseRegion>> exonsByChromosome = new HashMap<>();
        for(ContigEntry entry : entries)
        {
            List<BaseRegion> exons = exonsByChromosome.computeIfAbsent(entry.chromosome(), k -> new ArrayList<>());
            for(BaseRegion exon : entry.exonSpans())
                exons.add(new BaseRegion(exon.start(), exon.end()));
        }

        exonsByChromosome.values().forEach(ExonRegionIndex::mergeSpans);
        return new ExonRegionIndex(exonsByChromosome);
    }

    // Union-merge overlapping/adjacent spans into sorted, non-overlapping ranges. Deliberately not
    // BaseRegion.checkMergeOverlaps: that overwrites the end with the next span's (setEnd, not max), so an exon
    // nested inside an earlier-starting one shrinks the merged span - and exons nest across isoforms.
    private static void mergeSpans(final List<BaseRegion> regions)
    {
        regions.sort(Comparator.comparingInt(BaseRegion::start));

        int index = 0;
        while(index + 1 < regions.size())
        {
            BaseRegion current = regions.get(index);
            BaseRegion next = regions.get(index + 1);
            if(next.start() <= current.end() + 1)
            {
                current.setEnd(Math.max(current.end(), next.end()));
                regions.remove(index + 1);
            }
            else
            {
                ++index;
            }
        }
    }
}
