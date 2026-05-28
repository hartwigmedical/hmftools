package com.hartwig.hmftools.redux.splice;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Per-chromosome sorted exon spans, used by LiftBackResolver to decide whether a "hidden tie"
// (XS==AS with no XA) on a ref-only primary represents a sub-threshold paralog (rescue) or a real
// equally-scoring placement (keep MAPQ at 0). When the lifted primary's first base falls inside an
// annotated exon, the placement is likely the biologically-correct one and the tie's alt is almost
// certainly intronic/intergenic — rescue MAPQ. Outside-exon ties stay ambiguous.
//
// Loaded from the ensembl_data_cache CSVs (same files AnnotatedJunctionLoader reads). Lookup is
// O(log N) per chromosome via binary search over start positions.
public final class ExonRegionIndex
{
    private final Map<String, int[]> mExonStarts;
    private final Map<String, int[]> mExonEnds;

    private ExonRegionIndex(final Map<String, int[]> starts, final Map<String, int[]> ends)
    {
        mExonStarts = starts;
        mExonEnds = ends;
    }

    public boolean contains(final String chrom, final int pos)
    {
        final String key = normalize(chrom);
        final int[] starts = mExonStarts.get(key);
        final int[] ends = mExonEnds.get(key);
        if(starts == null)
            return false;

        // Intervals are merged at load time so they're sorted AND non-overlapping. A single binary
        // search for the largest start <= pos identifies the only possible containing range; the
        // earlier mistake was a backward walk that scanned every smaller-start interval since
        // starts[i] <= pos is trivially true for every i below the hit.
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

    public static ExonRegionIndex load(final String ensemblDir)
            throws IOException
    {
        final Map<String, String> geneChrom = loadGeneChromosomes(ensemblDir + "/ensembl_gene_data.csv");
        final Map<String, List<int[]>> byChrom = new HashMap<>();
        try(BufferedReader reader = new BufferedReader(new FileReader(ensemblDir + "/ensembl_trans_exon_data.csv")))
        {
            final String[] header = reader.readLine().split(",");
            final int gi = indexOf(header, "GeneId");
            final int si = indexOf(header, "ExonStart");
            final int ei = indexOf(header, "ExonEnd");
            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] cols = line.split(",");
                final String chrom = geneChrom.get(cols[gi]);
                if(chrom == null)
                    continue;
                byChrom.computeIfAbsent(chrom, k -> new ArrayList<>())
                        .add(new int[] { Integer.parseInt(cols[si]), Integer.parseInt(cols[ei]) });
            }
        }

        final Map<String, int[]> starts = new HashMap<>();
        final Map<String, int[]> ends = new HashMap<>();
        for(Map.Entry<String, List<int[]>> entry : byChrom.entrySet())
        {
            final List<int[]> merged = mergeIntervals(entry.getValue());
            final int[] s = new int[merged.size()];
            final int[] e = new int[merged.size()];
            for(int i = 0; i < merged.size(); ++i)
            {
                s[i] = merged.get(i)[0];
                e[i] = merged.get(i)[1];
            }
            starts.put(entry.getKey(), s);
            ends.put(entry.getKey(), e);
        }
        return new ExonRegionIndex(starts, ends);
    }

    // Sort by start, then sweep collapsing overlapping/abutting ranges. Result is sorted AND non-
    // overlapping, which is what makes the binary-search-only lookup safe (no backward walk needed).
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

    private static Map<String, String> loadGeneChromosomes(final String path)
            throws IOException
    {
        final Map<String, String> out = new HashMap<>();
        try(BufferedReader reader = new BufferedReader(new FileReader(path)))
        {
            final String[] header = reader.readLine().split(",");
            final int gi = indexOf(header, "GeneId");
            final int ci = indexOf(header, "Chromosome");
            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] cols = line.split(",");
                out.put(cols[gi], normalize(cols[ci]));
            }
        }
        return out;
    }

    // Ensembl CSVs ship bare chromosome names ("1") while V38 BAMs use "chr1". Normalize both
    // sides on the lookup so we don't care which convention the inputs use.
    private static String normalize(final String chrom)
    {
        if(chrom == null)
            return null;
        return chrom.startsWith("chr") ? chrom.substring(3) : chrom;
    }

    private static int indexOf(final String[] header, final String name)
    {
        for(int i = 0; i < header.length; ++i)
        {
            if(header[i].equals(name))
                return i;
        }
        throw new IllegalArgumentException("missing column: " + name);
    }
}
