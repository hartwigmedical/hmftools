package com.hartwig.hmftools.tars.liftback;

import com.hartwig.hmftools.tars.common.TarsConstants;

// One component of a record's alignment set: the record itself or one of its XA alts, post-lift.
// TransName/GeneId/GeneName are set only for Tx-contig alignments. SoftClipAtBoundary is true only
// when the leading/trailing S abuts an interior exon boundary (not the contig's outermost edge).
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
    // +1 forward / -1 reverse for tx-contig alignments; 0 otherwise (no transcript strand known).
    public final int TranscriptStrand;

    // Set by LiftBackResolver on the losing side of a discriminator decision. Kept for TSV-B diagnostics
    // but excluded from the BAM XA tag.
    public boolean Dropped = false;

    // Marks the alignment chosen as the BAM primary. Discriminator may flip it to an XA alt when bwa's
    // primary pick was suboptimal (e.g. processed pseudogene over spliced parent).
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

    // True iff at least one N is flanked by >= MIN_JUNCTION_ANCHOR M both sides (one strong junction trusts the read).
    public boolean cigarHasRealNJunction()
    {
        if(LiftedCigar == null)
        {
            return false;
        }

        int prevMLength = 0;
        int currentNum = 0;

        for(int i = 0; i < LiftedCigar.length(); ++i)
        {
            char c = LiftedCigar.charAt(i);
            if(c >= '0' && c <= '9')
            {
                currentNum = currentNum * 10 + (c - '0');
                continue;
            }

            if(c == 'N')
            {
                int nextM = nextMLength(LiftedCigar, i + 1);
                if(prevMLength >= TarsConstants.MIN_JUNCTION_ANCHOR && nextM >= TarsConstants.MIN_JUNCTION_ANCHOR)
                {
                    return true;
                }
                prevMLength = 0;
            }
            else if(c == 'M' || c == '=' || c == 'X')
            {
                prevMLength = currentNum;
            }
            else
            {
                // S/I/D/H/P break the adjacent anchor
                prevMLength = 0;
            }

            currentNum = 0;
        }

        return false;
    }

    private static int nextMLength(final String cigar, final int startIndex)
    {
        int num = 0;
        for(int i = startIndex; i < cigar.length(); ++i)
        {
            char c = cigar.charAt(i);
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
                return 0;
            }
        }
        return 0;
    }
}
