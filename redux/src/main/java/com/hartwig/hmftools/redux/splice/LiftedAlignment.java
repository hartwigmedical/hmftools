package com.hartwig.hmftools.redux.splice;

// one component of a record's alignment set: the record itself or one of its XA alts, post-lift.
// TransName / GeneId / GeneName are populated only when the source contig was a Tx contig.
// SoftClipAtBoundary is true only for Tx alignments whose leading or trailing S abuts an interior exon
// boundary in the source contig (an actual exon-junction edge, not the contig's outermost edge).
// ForwardStrand reflects the alignment's strand on the genome (Tx contigs are forward genomic, so strand
// passes through unchanged from the contig-space alignment).
public class LiftedAlignment
{
    public enum AlignmentSource
    {
        SELF,
        XA_INPUT
    }

    public final AlignmentSource Source;
    public final String OrigContig;
    public final int OrigPos;
    public final String OrigCigar;
    public final String LiftedChrom;
    public final int LiftedPos;
    public final String LiftedCigar;
    public final int AlignmentScore;
    public final int NumMismatches;
    public final String TransName;
    public final String GeneId;
    public final String GeneName;
    public final boolean SoftClipAtBoundary;
    public final boolean ForwardStrand;
    // Transcript's genome strand for tx-contig alignments (+1 forward / -1 reverse). 0 when this
    // alignment didn't come off a tx contig — the writer treats that as "no transcript strand known".
    public final int TranscriptStrand;

    // non-final: set by LiftBackResolver when this alignment is on the losing side of a confident
    // discriminator outcome (JUNCTION_FAVOURS_TX, TX_PICKED_REF_DROPPED, INTRON_RETAINED_FAVOURS_REF,
    // CROSS_LOCUS_FAVOURS_TX). Dropped alignments stay in LiftedAlignments for TSV-B diagnostics but are
    // excluded from the BAM XA tag.
    public boolean Dropped = false;

    // non-final: marks the alignment chosen as the BAM primary. Set true on SELF by default; the
    // discriminator may flip it to a winning XA alt when bwa's primary pick was suboptimal (e.g. processed
    // pseudogene preferred over the spliced parent gene). The displaced original becomes an XA entry.
    public boolean IsPrimaryChoice = false;

    public LiftedAlignment(
            final AlignmentSource source, final String origContig, final int origPos, final String origCigar,
            final String liftedChrom, final int liftedPos, final String liftedCigar,
            final int alignmentScore, final int numMismatches,
            final String transName, final String geneId, final String geneName,
            final boolean softClipAtBoundary, final boolean forwardStrand)
    {
        this(source, origContig, origPos, origCigar, liftedChrom, liftedPos, liftedCigar,
                alignmentScore, numMismatches, transName, geneId, geneName,
                softClipAtBoundary, forwardStrand, 0);
    }

    public LiftedAlignment(
            final AlignmentSource source, final String origContig, final int origPos, final String origCigar,
            final String liftedChrom, final int liftedPos, final String liftedCigar,
            final int alignmentScore, final int numMismatches,
            final String transName, final String geneId, final String geneName,
            final boolean softClipAtBoundary, final boolean forwardStrand, final int transcriptStrand)
    {
        Source = source;
        OrigContig = origContig;
        OrigPos = origPos;
        OrigCigar = origCigar;
        LiftedChrom = liftedChrom;
        LiftedPos = liftedPos;
        LiftedCigar = liftedCigar;
        AlignmentScore = alignmentScore;
        NumMismatches = numMismatches;
        TransName = transName;
        GeneId = geneId;
        GeneName = geneName;
        SoftClipAtBoundary = softClipAtBoundary;
        ForwardStrand = forwardStrand;
        TranscriptStrand = transcriptStrand;
    }

    public boolean fromTxContig()
    {
        return TransName != null;
    }

    public boolean cigarHasN()
    {
        return LiftedCigar != null && LiftedCigar.indexOf('N') >= 0;
    }

    public boolean cigarHasSoftClip()
    {
        return LiftedCigar != null && LiftedCigar.indexOf('S') >= 0;
    }

    // returns true iff the CIGAR has at least one N operator AND every N is flanked by an M
    // segment of at least minAnchor bp on both sides. Tiny anchors (e.g. 1-3M) on either side of
    // an N happen by chance when a tx-contig alignment's tail bases bleed across an exon boundary —
    // those are not real junctions and must not be used as evidence to swap the primary off a clean
    // ref full-match. The check walks the CIGAR once.
    public boolean cigarHasRealNJunction(final int minAnchor)
    {
        if(LiftedCigar == null)
            return false;

        boolean foundAny = false;
        int prevMLength = 0;
        int currentNum = 0;

        for(int i = 0; i < LiftedCigar.length(); ++i)
        {
            final char c = LiftedCigar.charAt(i);
            if(c >= '0' && c <= '9')
            {
                currentNum = currentNum * 10 + (c - '0');
                continue;
            }

            if(c == 'N')
            {
                if(prevMLength < minAnchor)
                    return false;
                // scan ahead for the next M-length; if no M follows, or it is below threshold, fail.
                final int nextM = nextMLength(LiftedCigar, i + 1);
                if(nextM < minAnchor)
                    return false;
                foundAny = true;
                prevMLength = 0;
            }
            else if(c == 'M' || c == '=' || c == 'X')
            {
                prevMLength = currentNum;
            }
            else
            {
                // any other op (S/I/D/H/P) breaks the "directly adjacent" anchor — treat as zero-length
                // M for the next N's purposes.
                prevMLength = 0;
            }

            currentNum = 0;
        }

        return foundAny;
    }

    private static int nextMLength(final String cigar, final int startIndex)
    {
        int num = 0;
        for(int i = startIndex; i < cigar.length(); ++i)
        {
            final char c = cigar.charAt(i);
            if(c >= '0' && c <= '9')
            {
                num = num * 10 + (c - '0');
            }
            else if(c == 'M' || c == '=' || c == 'X')
            {
                return num;
            }
            else
            {
                // any non-match op before we hit M means there is no directly adjacent M anchor.
                return 0;
            }
        }
        return 0;
    }
}
