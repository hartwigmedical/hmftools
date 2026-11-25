package com.hartwig.hmftools.tars.liftback;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.readToString;
import static com.hartwig.hmftools.tars.common.TarsCigarUtils.mergeAdjacentSameOp;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_EXTEND;
import static com.hartwig.hmftools.tars.common.TarsConstants.GAP_OPEN;
import static com.hartwig.hmftools.tars.common.TarsConstants.MATCH;
import static com.hartwig.hmftools.tars.common.TarsConstants.MISMATCH;

import static htsjdk.samtools.util.SequenceUtil.basesEqual;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import com.hartwig.hmftools.common.bam.CigarUtils;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.TextCigarCodec;

// Cigars are left unchanged: the only soft-clip walk in TARS is the weak-overhang collapse in OverhangGate.
public class GenomicAlignmentScorer
{
    private final RefGenomeInterface mRefGenome;

    public GenomicAlignmentScorer(final RefGenomeInterface refGenome)
    {
        mRefGenome = refGenome;
    }

    public void scorePlacements(final List<LiftedAlignment> alignments, final SAMRecord record)
    {
        if(alignments.size() < 2)
        {
            return;
        }

        ReadBases readBases = ReadBases.of(record);

        for(LiftedAlignment alignment : alignments)
        {
            if(alignment.Dropped)
            {
                continue;
            }

            alignment.GenomicScore = scoreAlignment(alignment, record, readBases);
        }
    }

    public int scoreRecord(final LiftedAlignment alignment, final SAMRecord record)
    {
        return scoreAlignment(alignment, record, ReadBases.of(record));
    }

    private int scoreAlignment(final LiftedAlignment alignment, final SAMRecord record, final ReadBases readBases)
    {
        byte[] alignmentBases = readBases.forAlignment(alignment);
        Cigar cigar = CigarUtils.cigarFromStr(alignment.LiftedCigar);

        int chromosomeLength = mRefGenome.getChromosomeLength(alignment.LiftedChromosome);
        int alignedEnd = alignment.LiftedPos + cigar.getReferenceLength() - 1;
        if(alignment.LiftedPos < 1 || alignedEnd > chromosomeLength)
        {
            throw new IllegalStateException(format("read(%s) alignment(%s:%d %s) outside chromosomeLength(%d)",
                    readToString(record), alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar,
                    chromosomeLength));
        }

        if(cigar.getReadLength() > alignmentBases.length)
        {
            throw new IllegalStateException(format("read(%s) alignment(%s:%d %s) cigarReadLength(%d) exceeds readLength(%d)",
                    readToString(record), alignment.LiftedChromosome, alignment.LiftedPos, alignment.LiftedCigar,
                    cigar.getReadLength(), alignmentBases.length));
        }

        return genomicAlignmentScore(mRefGenome, alignment.LiftedChromosome, alignment.LiftedPos, cigar, alignmentBases);
    }

    // An intron-only lift preserves BWA's scoring operations: the transcript contig contains the same exonic bases and
    // N has no alignment-score cost. Any other CIGAR change makes the input AS stale and requires a genomic rescore.
    public static boolean requiresRescore(final SAMRecord record, final LiftedAlignment alignment)
    {
        List<CigarElement> liftedOps = new ArrayList<>();
        for(CigarElement element : TextCigarCodec.decode(alignment.LiftedCigar).getCigarElements())
        {
            if(element.getOperator() != CigarOperator.N)
            {
                liftedOps.add(element);
            }
        }

        List<CigarElement> inputOps = new ArrayList<>(record.getCigar().getCigarElements());
        if(record.getReadNegativeStrandFlag() == alignment.ForwardStrand)
        {
            Collections.reverse(inputOps);
        }

        return !mergeAdjacentSameOp(liftedOps).equals(mergeAdjacentSameOp(inputOps));
    }

    public static boolean requiresRescore(final SAMRecord record, final LiftedRecord lifted)
    {
        return lifted.primaryIndex() != 0 || requiresRescore(record, lifted.primaryAlignment());
    }

    // Null unless the whole window is present, so nothing is scored against truncated bases: past a contig end
    // RefGenomeSource lets htsjdk throw while MockRefGenome returns empty.
    @Nullable
    public static byte[] referenceBases(
            final RefGenomeInterface refGenome, final String chromosome, final int posStart, final int posEnd)
    {
        if(posStart < 1 || posEnd < posStart || posEnd > refGenome.getChromosomeLength(chromosome))
        {
            return null;
        }

        byte[] bases = refGenome.getBases(chromosome, posStart, posEnd);
        return bases != null && bases.length == posEnd - posStart + 1 ? bases : null;
    }

    public static int alignedBaseScore(final byte readBase, final byte refBase)
    {
        return basesEqual(readBase, refBase) ? MATCH : MISMATCH;
    }

    // bwa-mem score of read vs ref over their overlapping length.
    public static int ungappedScore(final byte[] read, final byte[] ref)
    {
        int total = 0;
        int length = Math.min(read.length, ref.length);
        for(int i = 0; i < length; ++i)
        {
            total += alignedBaseScore(read[i], ref[i]);
        }
        return total;
    }

    // Longest prefix reaching the highest cumulative bwa-mem score, 0 if none positive; the >= tie extends past an internal mismatch.
    public static int bestScoringPrefixLength(final byte[] read, final byte[] ref)
    {
        int length = Math.min(read.length, ref.length);
        int score = 0;
        int bestScore = 0;
        int bestLength = 0;
        for(int i = 0; i < length; ++i)
        {
            score += alignedBaseScore(read[i], ref[i]);
            if(score >= bestScore && score > 0)
            {
                bestScore = score;
                bestLength = i + 1;
            }
        }
        return bestLength;
    }

    public static int genomicAlignmentScore(
            final RefGenomeInterface refGenome, final String chromosome, final int alignmentStart, final Cigar cigar,
            final byte[] readBases)
    {
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
                    byte[] ref = referenceBases(refGenome, chromosome, refPos, refPos + length - 1);
                    for(int i = 0; i < length; ++i)
                    {
                        totalScore += alignedBaseScore(readBases[queryPos + i], ref[i]);
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
