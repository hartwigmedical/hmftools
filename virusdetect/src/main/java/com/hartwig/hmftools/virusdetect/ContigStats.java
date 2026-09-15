package com.hartwig.hmftools.virusdetect;

import java.util.Optional;

// Per-contig support over a set of aligned reads, plus how it fares against the sibling contigs of the same virus.
public record ContigStats(
        String contig,
        int contigLength,
        // reads with any alignment to this contig
        int readCount,
        // reads with more than one alignment to this contig (BWA -a repeats/multi-loci)
        int multiAlignReads,
        // alignments to this contig per read (>= 1)
        SummaryStats alignPerRead,
        // alignments dropped for clipping over the contig end (circular-genome artifact)
        int originClippedReads,
        // contig positions with at least one aligned base
        int coveredBases,
        // depth spans the whole contig, so uncovered positions count as depth 0
        SummaryStats depth,
        // BWA alignment score distribution across those reads
        SummaryStats alignerScore,
        // strain rivalry: reads softly attributed to this strain, split across the sibling contigs by divergence
        double readVotes,
        // reads that strictly beat every rival contig (ties credit no contig)
        int readsBestInRivals,
        // divergence lead over the runner-up on strictly-won reads; empty when there were none
        Optional<SummaryStats> margins)
{
    public double coverageFraction()
    {
        return contigLength == 0 ? 0 : (double) coveredBases / contigLength;
    }
}
