package com.hartwig.hmftools.tars.liftback;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.ExonData;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

// Per-chromosome merged exon spans loaded from the ensembl data cache. Used by LiftBackResolver to
// distinguish a hidden tie where the primary lands in an annotated exon (rescue MAPQ) from a genuinely
// ambiguous tie (keep MAPQ 0). Lookup is O(log N) via binary search over sorted, non-overlapping spans.
// Chromosomes are keyed in the run's ref-genome form (versionedChromosome) so lookups by lifted genomic
// contig match directly.
public final class ExonRegionIndex
{
    private final Map<String, int[]> mExonStarts;
    private final Map<String, int[]> mExonEnds;

    private ExonRegionIndex(final Map<String, int[]> starts, final Map<String, int[]> ends)
    {
        mExonStarts = starts;
        mExonEnds = ends;
    }

    public boolean contains(final String chromosome, final int pos)
    {
        final int[] starts = mExonStarts.get(chromosome);
        final int[] ends = mExonEnds.get(chromosome);
        if(starts == null)
            return false;

        // Binary search for the largest start <= pos; merged intervals guarantee only one candidate.
        int lo = 0;
        int hi = starts.length - 1;
        int idx = -1;
        while(lo <= hi)
        {
            final int mid = (lo + hi) >>> 1;
            if(starts[mid] <= pos)
            {
                idx = mid;
                lo = mid + 1;
            }
            else
            {
                hi = mid - 1;
            }
        }
        return idx >= 0 && pos <= ends[idx];
    }

    public static ExonRegionIndex load(final String ensemblDir, final RefGenomeVersion refGenomeVersion)
    {
        final EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, false);
        ensemblDataCache.load(false);
        return fromCache(ensemblDataCache, refGenomeVersion);
    }

    public static ExonRegionIndex fromCache(final EnsemblDataCache ensemblDataCache, final RefGenomeVersion refGenomeVersion)
    {
        final Map<String, List<int[]>> spansByChromosome = new HashMap<>();
        for(final List<GeneData> genes : ensemblDataCache.getChrGeneDataMap().values())
        {
            for(final GeneData gene : genes)
            {
                final List<TranscriptData> transcripts = ensemblDataCache.getTranscripts(gene.GeneId);
                if(transcripts == null)
                    continue;

                final String chromosome = refGenomeVersion.versionedChromosome(gene.Chromosome);
                final List<int[]> spans = spansByChromosome.computeIfAbsent(chromosome, k -> new ArrayList<>());
                for(final TranscriptData transcript : transcripts)
                {
                    for(final ExonData exon : transcript.exons())
                        spans.add(new int[] { exon.Start, exon.End });
                }
            }
        }

        final Map<String, int[]> starts = new HashMap<>();
        final Map<String, int[]> ends = new HashMap<>();
        for(final Map.Entry<String, List<int[]>> entry : spansByChromosome.entrySet())
        {
            final List<int[]> merged = mergeIntervals(entry.getValue());
            final int[] mergedStarts = new int[merged.size()];
            final int[] mergedEnds = new int[merged.size()];
            for(int i = 0; i < merged.size(); ++i)
            {
                mergedStarts[i] = merged.get(i)[0];
                mergedEnds[i] = merged.get(i)[1];
            }
            starts.put(entry.getKey(), mergedStarts);
            ends.put(entry.getKey(), mergedEnds);
        }
        return new ExonRegionIndex(starts, ends);
    }

    // Produces sorted, non-overlapping intervals required for safe binary-search lookup.
    private static List<int[]> mergeIntervals(final List<int[]> intervals)
    {
        if(intervals.isEmpty())
            return intervals;

        Collections.sort(intervals, (a, b) -> Integer.compare(a[0], b[0]));
        final List<int[]> merged = new ArrayList<>();
        int[] current = new int[] { intervals.get(0)[0], intervals.get(0)[1] };
        for(int i = 1; i < intervals.size(); ++i)
        {
            final int[] next = intervals.get(i);
            if(next[0] <= current[1] + 1)
                current[1] = Math.max(current[1], next[1]);
            else
            {
                merged.add(current);
                current = new int[] { next[0], next[1] };
            }
        }
        merged.add(current);
        return merged;
    }
}
