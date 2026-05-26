package com.hartwig.hmftools.redux.splice.rescue;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

// Loads annotated splice junctions from the ensembl_data_cache CSVs. Same source the python
// BamCompare analysis script reads. Returns a Set<ChrIntron> keyed by (chrom, intronStart,
// intronEnd) so the rescue resolver can do O(1) annotated-junction lookups.
//
// Algorithm: pair each transcript's adjacent exons in TransId+Rank order and emit the intron
// between them. Multiple transcripts may share the same junction; the Set deduplicates naturally.
public final class AnnotatedJunctionLoader
{
    private AnnotatedJunctionLoader() {}

    public static Set<ChrIntron> load(final String ensemblDir)
            throws IOException
    {
        final Map<String, String> geneChrom = loadGeneChromosomes(ensemblDir + "/ensembl_gene_data.csv");
        final Map<String, List<Exon>> exonsByTrans = loadExonsByTrans(ensemblDir + "/ensembl_trans_exon_data.csv");

        final Set<ChrIntron> introns = new HashSet<>();
        for(Map.Entry<String, List<Exon>> entry : exonsByTrans.entrySet())
        {
            final List<Exon> exons = entry.getValue();
            if(exons.size() < 2)
                continue;

            Collections.sort(exons, Comparator.comparingInt(e -> e.Rank));
            final String chrom = geneChrom.get(exons.get(0).GeneId);
            if(chrom == null)
                continue;

            for(int i = 0; i < exons.size() - 1; ++i)
            {
                final Exon a = exons.get(i);
                final Exon b = exons.get(i + 1);
                // exons may be ranked in transcription order regardless of genome strand. The
                // intron lies between the higher genomic end and the lower genomic start of
                // adjacent exons — work it out from coords directly.
                final int prevEnd = Math.min(a.End, b.End);
                final int nextStart = Math.max(a.Start, b.Start);
                if(nextStart <= prevEnd + 1)
                    continue;     // no gap → not an intron
                final int intronStart = Math.min(a.End, b.End) + 1;
                final int intronEnd = Math.max(a.Start, b.Start) - 1;
                if(intronEnd < intronStart)
                    continue;
                introns.add(new ChrIntron(chrom, intronStart, intronEnd));
            }
        }

        return introns;
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
                final String chrom = cols[ci].startsWith("chr") ? cols[ci] : "chr" + cols[ci];
                out.put(cols[gi], chrom);
            }
        }
        return out;
    }

    private static Map<String, List<Exon>> loadExonsByTrans(final String path)
            throws IOException
    {
        final Map<String, List<Exon>> out = new HashMap<>();
        try(BufferedReader reader = new BufferedReader(new FileReader(path)))
        {
            final String[] header = reader.readLine().split(",");
            final int gi = indexOf(header, "GeneId");
            final int ti = indexOf(header, "TransId");
            final int si = indexOf(header, "ExonStart");
            final int ei = indexOf(header, "ExonEnd");
            final int ri = indexOf(header, "ExonRank");
            String line;
            while((line = reader.readLine()) != null)
            {
                final String[] cols = line.split(",");
                final Exon exon = new Exon(
                        cols[gi], Integer.parseInt(cols[si]), Integer.parseInt(cols[ei]),
                        Integer.parseInt(cols[ri]));
                out.computeIfAbsent(cols[ti], k -> new ArrayList<>()).add(exon);
            }
        }
        return out;
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

    private static final class Exon
    {
        final String GeneId;
        final int Start;
        final int End;
        final int Rank;

        Exon(final String geneId, final int start, final int end, final int rank)
        {
            GeneId = geneId;
            Start = start;
            End = end;
            Rank = rank;
        }
    }
}
