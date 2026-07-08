package com.hartwig.hmftools.tars.liftback;

import com.hartwig.hmftools.tars.common.TarsConstants;
import com.hartwig.hmftools.tars.liftback.supplementary.RefSequenceSource;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.TextCigarCodec;

// bwa-mem-style score of a lifted alignment against the reference genome: match/mismatch over aligned bases, an
// affine gap for each I/D, introns (N) and clips (S/H) free. bwa's own AlignmentScore is scored in tx-contig space
// for a transcript candidate and is absent for XA alts, so it cannot rank a read's candidate placements against
// each other; this recomputes every candidate on one comparable genomic scale for the Step 2 tie-break.
public final class AlignmentScorer
{
    private AlignmentScorer() { }

    // readBases must already be in the alignment's own orientation (reverse-complemented for an opposite-strand alt).
    // Returns Integer.MIN_VALUE when the reference is unavailable, matching LiftedAlignment's "not computed" sentinel.
    public static int score(
            final RefSequenceSource refSource, final String chromosome, final int alignmentStart,
            final String cigarStr, final byte[] readBases)
    {
        if(refSource == null || readBases == null || cigarStr == null)
        {
            return Integer.MIN_VALUE;
        }

        Cigar cigar = TextCigarCodec.decode(cigarStr);
        int score = 0;
        int queryPos = 0;
        int refPos = alignmentStart;

        for(final CigarElement element : cigar.getCigarElements())
        {
            int length = element.getLength();
            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    byte[] ref = refSource.getBases(chromosome, refPos, refPos + length - 1);
                    for(int i = 0; i < length; ++i)
                    {
                        boolean match = ref != null && i < ref.length && queryPos + i < readBases.length
                                && basesEqualIgnoreCase(readBases[queryPos + i], ref[i]);
                        score += match ? TarsConstants.MATCH : TarsConstants.MISMATCH;
                    }
                    queryPos += length;
                    refPos += length;
                    break;

                case I:
                    score += TarsConstants.GAP_OPEN + (length - 1) * TarsConstants.GAP_EXTEND;
                    queryPos += length;
                    break;

                case D:
                    score += TarsConstants.GAP_OPEN + (length - 1) * TarsConstants.GAP_EXTEND;
                    refPos += length;
                    break;

                case N:
                    refPos += length;
                    break;

                case S:
                    queryPos += length;
                    break;

                default: // H, P consume neither aligned query nor scored reference
                    break;
            }
        }

        return score;
    }

    private static boolean basesEqualIgnoreCase(final byte a, final byte b)
    {
        return (a & ~0x20) == (b & ~0x20);
    }
}
