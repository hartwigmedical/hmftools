package com.hartwig.hmftools.tars.common;

import static com.hartwig.hmftools.tars.common.TarsConstants.MATCH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MISMATCH;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_EXTEND;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_OPEN;

import static htsjdk.samtools.util.SequenceUtil.basesEqual;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;

// Scores read bases against the reference the way bwa-mem does, for the passes that re-check a lifted placement.
public final class BwaScoring
{
    private BwaScoring() { }

    // Null unless the whole requested window is present, so nothing is ever scored against truncated bases. The
    // implementations disagree here: RefGenomeSource lets htsjdk throw past a contig end, MockRefGenome returns empty.
    public static byte[] refWindow(
            final RefGenomeInterface refGenome, final String chromosome, final int posStart, final int posEnd)
    {
        if(refGenome == null || posStart < 1 || posEnd < posStart || posEnd > refGenome.getChromosomeLength(chromosome))
        {
            return null;
        }

        byte[] bases = refGenome.getBases(chromosome, posStart, posEnd);
        return bases != null && bases.length == posEnd - posStart + 1 ? bases : null;
    }

    // one base pair's contribution to a bwa-mem score
    public static int baseScore(final byte readBase, final byte refBase)
    {
        return basesEqual(readBase, refBase) ? MATCH : MISMATCH;
    }

    // Straight bwa-mem score of read vs ref over their overlapping length.
    public static int score(final byte[] read, final byte[] ref)
    {
        int total = 0;
        int length = Math.min(read.length, ref.length);
        for(int i = 0; i < length; ++i)
        {
            total += baseScore(read[i], ref[i]);
        }
        return total;
    }

    // Longest prefix reaching the highest cumulative bwa-mem score (0 if none positive); the >= tie extends past an
    // internal mismatch.
    public static int maxScoringPrefix(final byte[] read, final byte[] ref)
    {
        int length = Math.min(read.length, ref.length);
        int score = 0;
        int bestScore = 0;
        int bestLength = 0;
        for(int i = 0; i < length; ++i)
        {
            score += baseScore(read[i], ref[i]);
            if(score >= bestScore && score > 0)
            {
                bestScore = score;
                bestLength = i + 1;
            }
        }
        return bestLength;
    }

    public static int genomicScore(
            final RefGenomeInterface refGenome, final String chromosome, final int alignmentStart, final String cigarStr,
            final byte[] readBases)
    {
        if(refGenome == null || readBases == null || cigarStr == null)
        {
            return Integer.MIN_VALUE;
        }

        Cigar cigar = CigarUtils.cigarFromStr(cigarStr);
        int totalScore = 0;
        int queryPos = 0;
        int refPos = alignmentStart;

        for(CigarElement element : cigar.getCigarElements())
        {
            int length = element.getLength();
            switch(element.getOperator())
            {
                case M:
                case EQ:
                case X:
                    byte[] ref = refWindow(refGenome, chromosome, refPos, refPos + length - 1);
                    for(int i = 0; i < length; ++i)
                    {
                        boolean scorable = ref != null && queryPos + i < readBases.length;
                        totalScore += scorable ? baseScore(readBases[queryPos + i], ref[i]) : MISMATCH;
                    }
                    queryPos += length;
                    refPos += length;
                    break;

                case I:
                    totalScore += GAP_OPEN + (length - 1) * GAP_EXTEND;
                    queryPos += length;
                    break;

                case D:
                    totalScore += GAP_OPEN + (length - 1) * GAP_EXTEND;
                    refPos += length;
                    break;

                case N:
                    refPos += length;
                    break;

                case S:
                    queryPos += length;
                    break;

                default:
                    break;
            }
        }

        return totalScore;
    }
}
