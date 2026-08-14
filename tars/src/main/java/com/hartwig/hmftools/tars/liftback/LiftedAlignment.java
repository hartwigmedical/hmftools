package com.hartwig.hmftools.tars.liftback;

import static com.hartwig.hmftools.common.region.BaseRegion.positionsOverlap;

import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.region.ChrBaseRegion;

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
    // +1 forward / -1 reverse for tx-contig alignments; 0 when no transcript strand is known.
    public final int TranscriptStrand;
    public final List<Integer> MergedSupplementaryIndices;
    public final List<ChrBaseRegion> MergedSupplementaryIntrons;
    public final int MergedSupplementaryMapQuality;

    // set when an overhang collapse made this alt a duplicate placement; excluded from the XA tag.
    public boolean Dropped = false;

    // bwa-mem-style score of the lifted placement against the genome, letting tx and ref candidates compare on one scale.
    // Integer.MIN_VALUE = not computed.
    public int GenomicScore = Integer.MIN_VALUE;

    public LiftedAlignment(
            final String liftedChrom, final int liftedPos, final String liftedCigar, final int numMismatches,
            final boolean fromTxContig, final boolean softClipAtBoundary, final boolean forwardStrand, final int transcriptStrand)
    {
        this(
                liftedChrom, liftedPos, liftedCigar, numMismatches, fromTxContig, softClipAtBoundary, forwardStrand,
                transcriptStrand, Collections.emptyList(), Collections.emptyList(), -1);
    }

    private LiftedAlignment(
            final String liftedChrom, final int liftedPos, final String liftedCigar, final int numMismatches,
            final boolean fromTxContig, final boolean softClipAtBoundary, final boolean forwardStrand, final int transcriptStrand,
            final List<Integer> mergedSupplementaryIndices, final List<ChrBaseRegion> mergedSupplementaryIntrons,
            final int mergedSupplementaryMapQuality)
    {
        LiftedChromosome = liftedChrom;
        LiftedPos = liftedPos;
        LiftedCigar = liftedCigar;
        NumMismatches = numMismatches;
        FromTxContig = fromTxContig;
        SoftClipAtBoundary = softClipAtBoundary;
        ForwardStrand = forwardStrand;
        TranscriptStrand = transcriptStrand;
        MergedSupplementaryIndices = List.copyOf(mergedSupplementaryIndices);
        MergedSupplementaryIntrons = List.copyOf(mergedSupplementaryIntrons);
        MergedSupplementaryMapQuality = mergedSupplementaryMapQuality;
    }

    // copy with a revised lifted position and cigar; the boundary-softclip flag drops if the new cigar has no softclip.
    // the chosen primary is tracked by index on LiftedRecord, so no copy here carries it.
    public LiftedAlignment withLiftedCigar(final int liftedPos, final String liftedCigar)
    {
        boolean stillSoftClipped = liftedCigar.indexOf('S') >= 0;
        LiftedAlignment revised = new LiftedAlignment(
                LiftedChromosome, liftedPos, liftedCigar, NumMismatches,
                FromTxContig, SoftClipAtBoundary && stillSoftClipped, ForwardStrand, TranscriptStrand,
                MergedSupplementaryIndices, MergedSupplementaryIntrons, MergedSupplementaryMapQuality);
        revised.Dropped = Dropped;
        revised.GenomicScore = GenomicScore;
        return revised;
    }

    public LiftedAlignment withSupplementaryMerge(
            final int liftedPos, final String liftedCigar, final List<Integer> mergedSupplementaryIndices,
            final List<ChrBaseRegion> mergedSupplementaryIntrons, final int mergedSupplementaryMapQuality,
            final int spliceStrand)
    {
        int mergedTranscriptStrand = TranscriptStrand != 0 ? TranscriptStrand : spliceStrand;
        LiftedAlignment revised = new LiftedAlignment(
                LiftedChromosome, liftedPos, liftedCigar, NumMismatches,
                FromTxContig, false, ForwardStrand, mergedTranscriptStrand,
                mergedSupplementaryIndices, mergedSupplementaryIntrons, mergedSupplementaryMapQuality);
        revised.GenomicScore = GenomicScore;
        return revised;
    }

    public LiftedAlignment withTranscriptStrand(final int transcriptStrand)
    {
        LiftedAlignment revised = new LiftedAlignment(
                LiftedChromosome, LiftedPos, LiftedCigar, NumMismatches,
                FromTxContig, SoftClipAtBoundary, ForwardStrand, transcriptStrand,
                MergedSupplementaryIndices, MergedSupplementaryIntrons, MergedSupplementaryMapQuality);
        revised.Dropped = Dropped;
        revised.GenomicScore = GenomicScore;
        return revised;
    }

    public boolean hasSupplementaryMerge()
    {
        return !MergedSupplementaryIndices.isEmpty();
    }

    // bwa XA entry "chrom,<sign>pos,cigar,NM;", the position sign carrying the strand.
    public String toXaEntry()
    {
        return LiftedChromosome + ',' + (ForwardStrand ? '+' : '-') + LiftedPos + ','
                + LiftedCigar + ',' + NumMismatches + ';';
    }

    public String locusKey()
    {
        return LiftedChromosome + ":" + LiftedPos;
    }

    public AlignmentKey key()
    {
        return AlignmentKey.from(this);
    }

    public int alignedEnd()
    {
        return LiftedPos + CigarUtils.calcCigarAlignedLength(LiftedCigar) - 1;
    }

    // Same-locus test shared by the NH locus count and the XA build: the two must agree, or a read is emitted with
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
