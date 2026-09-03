package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import com.hartwig.hmftools.common.bam.CigarUtils;

// One component of a record's alignment set (the record itself or an XA alt), post-lift. SoftClipAtBoundary is set only
// when a leading/trailing S abuts an interior (not outermost) exon boundary.
public class LiftedAlignment
{
    public final String LiftedChromosome;
    public final int LiftedPos;
    public final String LiftedCigar;
    public final int NumMismatches;
    public final boolean FromTxContig;
    public final boolean SoftClipAtBoundary;
    public final boolean ForwardStrand;
    // +1 forward / -1 reverse for tx-contig alignments; 0 otherwise (no transcript strand known).
    public final int TranscriptStrand;

    // set when an overhang collapse made this alt a duplicate placement; excluded from the XA tag.
    public boolean Dropped = false;

    // bwa-mem-style score of the lifted placement against the genome, letting tx and ref candidates compare on one scale.
    // Integer.MIN_VALUE = not computed.
    public int GenomicScore = Integer.MIN_VALUE;

    public LiftedAlignment(
            final String liftedChrom, final int liftedPos, final String liftedCigar, final int numMismatches,
            final boolean fromTxContig, final boolean softClipAtBoundary, final boolean forwardStrand, final int transcriptStrand)
    {
        LiftedChromosome = liftedChrom;
        LiftedPos = liftedPos;
        LiftedCigar = liftedCigar;
        NumMismatches = numMismatches;
        FromTxContig = fromTxContig;
        SoftClipAtBoundary = softClipAtBoundary;
        ForwardStrand = forwardStrand;
        TranscriptStrand = transcriptStrand;
    }

    // copy with a revised lifted position and cigar; the boundary-softclip flag drops if the new cigar has no softclip.
    // Dropped and GenomicScore carry over; the chosen primary is tracked by index on LiftedRecord, not here.
    public LiftedAlignment withLiftedCigar(final int liftedPos, final String liftedCigar)
    {
        boolean stillSoftClipped = liftedCigar.indexOf('S') >= 0;
        LiftedAlignment revised = new LiftedAlignment(
                LiftedChromosome, liftedPos, liftedCigar, NumMismatches,
                FromTxContig, SoftClipAtBoundary && stillSoftClipped, ForwardStrand, TranscriptStrand);
        revised.Dropped = Dropped;
        revised.GenomicScore = GenomicScore;
        return revised;
    }

    // this placement as a bwa XA entry: "chrom,<sign>pos,cigar,NM;", the position sign carrying the strand.
    public String toXaEntry()
    {
        return LiftedChromosome + ',' + (ForwardStrand ? '+' : '-') + LiftedPos + ','
                + LiftedCigar + ',' + NumMismatches + ';';
    }

    // chrom:pos key grouping alignments by genomic locus.
    public String locusKey()
    {
        return LiftedChromosome + ":" + LiftedPos;
    }

    public int alignedEnd()
    {
        return LiftedPos + CigarUtils.calcCigarAlignedLength(LiftedCigar) - 1;
    }

    // Same-locus test shared by the NH locus count and the XA build - the two must agree, or a read is emitted with
    // NH=1 alongside a non-empty XA. False against a null other, so an absent primary collapses nothing.
    public boolean overlaps(final LiftedAlignment other)
    {
        if(other == null || !LiftedChromosome.equals(other.LiftedChromosome))
        {
            return false;
        }
        return positionsOverlap(LiftedPos, alignedEnd(), other.LiftedPos, other.alignedEnd());
    }

    public boolean cigarHasN()
    {
        return LiftedCigar != null && LiftedCigar.indexOf('N') >= 0;
    }

    public boolean cigarHasSoftClip()
    {
        return LiftedCigar != null && LiftedCigar.indexOf('S') >= 0;
    }
}
