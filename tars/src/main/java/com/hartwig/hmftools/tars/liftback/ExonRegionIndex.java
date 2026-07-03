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
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.tars.common.ContigEntry;

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
        int[] starts = mExonStarts.get(chromosome);
        int[] ends = mExonEnds.get(chromosome);
        if(starts == null)
        {
            return false;
        }

        // Binary search for the largest start <= pos; merged intervals guarantee only one candidate.
        int lo = 0;
        int hi = starts.length - 1;
        int idx = -1;
        while(lo <= hi)
        {
            int mid = (lo + hi) >>> 1;
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
        EnsemblDataCache ensemblDataCache = new EnsemblDataCache(ensemblDir, refGenomeVersion);
        ensemblDataCache.setRequiredData(true, false, false, false);
        ensemblDataCache.load(false);
        return fromCache(ensemblDataCache, refGenomeVersion);
    }

    // built from the contig sidecar instead of the ensembl cache: each entry's exonSpans are the genomic exons
    // of a multi-exon transcript (the only transcripts that get a tx-contig). Single-exon transcripts are absent,
    // but their exons only matter to the discriminator when a read has a tx-contig match -- which only multi-exon
    // transcripts produce -- so the derived index is behaviour-equivalent to the ensembl-loaded one. Keyed by the
    // sidecar chromosome, the exact form the lift emits.
    public static ExonRegionIndex fromContigEntries(final List<ContigEntry> entries)
    {
        Map<String, List<int[]>> spansByChromosome = new HashMap<>();
        for(final ContigEntry entry : entries)
        {
            List<int[]> spans = spansByChromosome.computeIfAbsent(entry.chromosome(), k -> new ArrayList<>());
            for(final BaseRegion exon : entry.exonSpans())
            {
                spans.add(new int[] { exon.start(), exon.end() });
            }
        }
        return fromSpans(spansByChromosome);
    }

    public static ExonRegionIndex fromCache(final EnsemblDataCache ensemblDataCache, final RefGenomeVersion refGenomeVersion)
    {
        Map<String, List<int[]>> spansByChromosome = new HashMap<>();
        for(final List<GeneData> genes : ensemblDataCache.getChrGeneDataMap().values())
        {
            for(final GeneData gene : genes)
            {
                List<TranscriptData> transcripts = ensemblDataCache.getTranscripts(gene.GeneId);
                if(transcripts == null)
                    continue;

                String chromosome = refGenomeVersion.versionedChromosome(gene.Chromosome);
                List<int[]> spans = spansByChromosome.computeIfAbsent(chromosome, k -> new ArrayList<>());
                for(final TranscriptData transcript : transcripts)
                {
                    for(final ExonData exon : transcript.exons())
                    {
                        spans.add(new int[] { exon.Start, exon.End });
                    }
                }
            }
        }
        return fromSpans(spansByChromosome);
    }

    private static ExonRegionIndex fromSpans(final Map<String, List<int[]>> spansByChromosome)
    {
        Map<String, int[]> starts = new HashMap<>();
        Map<String, int[]> ends = new HashMap<>();
        for(final Map.Entry<String, List<int[]>> entry : spansByChromosome.entrySet())
        {
            List<int[]> merged = mergeIntervals(entry.getValue());
            int[] mergedStarts = new int[merged.size()];
            int[] mergedEnds = new int[merged.size()];
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
        {
            return intervals;
        }

        Collections.sort(intervals, (a, b) -> Integer.compare(a[0], b[0]));
        List<int[]> merged = new ArrayList<>();
        int[] current = new int[] { intervals.get(0)[0], intervals.get(0)[1] };
        for(int i = 1; i < intervals.size(); ++i)
        {
            int[] next = intervals.get(i);
            if(next[0] <= current[1] + 1)
            {
                current[1] = Math.max(current[1], next[1]);
            }
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
