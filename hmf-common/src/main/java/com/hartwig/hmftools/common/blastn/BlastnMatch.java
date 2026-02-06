package com.hartwig.hmftools.common.blastn;

import static java.lang.Character.isDigit;
import static java.lang.Character.isLetter;
import static java.lang.Integer.parseInt;
import static java.lang.String.format;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;
import java.util.Optional;
import java.util.function.BiConsumer;

import com.hartwig.hmftools.common.genome.region.Strand;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;

// Represents a BLASTn match with various alignment details.
public class BlastnMatch
{
    public final int QuerySeqLen;
    public final String SubjectTitle;
    public final double PercentageIdent;
    public final double QueryCoverage;
    public final int AlignmentLength;
    public final int NumMismatch;
    public final int NumGapOpenings;
    public final int SubjectAlignStart;
    public final int SubjectAlignEnd;
    public final Strand SubjectFrame;
    public final double ExpectedValue;
    public final double BitScore;
    public final String AlignedPartOfQuerySeq;
    public final String AlignedPartOfSubjectSeq;
    public final String Operations;

    private int mQueryAlignStart;
    private int mQueryAlignEnd;

    public BlastnMatch(
            int querySeqLen, String subjectTitle, double percentageIdent, double queryCoverage, int alignmentLength,
            int numMismatch, int numGapOpenings, int queryAlignStart, int queryAlignEnd, int subjectAlignStart, int subjectAlignEnd,
            Strand subjectFrame, double expectedValue, double bitScore, String alignedPartOfQuerySeq, String alignedPartOfSubjectSeq,
            String operations)
    {
        QuerySeqLen = querySeqLen;
        SubjectTitle = subjectTitle;
        PercentageIdent = percentageIdent;
        QueryCoverage = queryCoverage;
        AlignmentLength = alignmentLength;
        NumMismatch = numMismatch;
        NumGapOpenings = numGapOpenings;
        mQueryAlignStart = queryAlignStart;
        mQueryAlignEnd = queryAlignEnd;
        SubjectAlignStart = subjectAlignStart;
        SubjectAlignEnd = subjectAlignEnd;
        SubjectFrame = subjectFrame;
        ExpectedValue = expectedValue;
        BitScore = bitScore;
        AlignedPartOfQuerySeq = alignedPartOfQuerySeq;
        AlignedPartOfSubjectSeq = alignedPartOfSubjectSeq;
        Operations = operations;
    }

    public int getQuerySeqLen()
    {
        return QuerySeqLen;
    }

    public String getSubjectTitle()
    {
        return SubjectTitle;
    }

    public double getPercentageIdent()
    {
        return PercentageIdent;
    }

    public double getQueryCoverage()
    {
        return QueryCoverage;
    }

    public int getAlignmentLength()
    {
        return AlignmentLength;
    }

    public int getNumMismatch()
    {
        return NumMismatch;
    }

    public int getNumGapOpenings()
    {
        return NumGapOpenings;
    }

    public int getQueryAlignStart()
    {
        return mQueryAlignStart;
    }
    public void setQueryAlignStart(final int queryAlignStart)
    {
        mQueryAlignStart = queryAlignStart;
    }

    public int getQueryAlignEnd()
    {
        return mQueryAlignEnd;
    }
    public void setQueryAlignEnd(final int queryAlignEnd) { mQueryAlignEnd = queryAlignEnd; }

    public int getSubjectAlignStart()
    {
        return SubjectAlignStart;
    }

    public int getSubjectAlignEnd()
    {
        return SubjectAlignEnd;
    }

    public Strand getSubjectFrame()
    {
        return SubjectFrame;
    }

    public double getExpectedValue()
    {
        return ExpectedValue;
    }

    public double getBitScore()
    {
        return BitScore;
    }

    public String getAlignedPartOfQuerySeq()
    {
        return AlignedPartOfQuerySeq;
    }

    public String getAlignedPartOfSubjectSeq()
    {
        return AlignedPartOfSubjectSeq;
    }

    public List<CigarElement> getCigar()
    {
        // Example operations format:
        // 10AG5T-3-C8
        // 10 matched, 1 mismatch A->G, 5 matched, 1 inserted T, 3 matched, 1 deleted C, 8 matched.

        ArrayList<CigarElement> cigar = new ArrayList<>();

        BiConsumer<CigarOperator, Integer> addOperation = (operator, length) -> {
            CigarElement prev = cigar.isEmpty() ? null : cigar.get(cigar.size() - 1);
            if (prev != null && prev.getOperator() == operator)
            {
                cigar.set(cigar.size() - 1, new CigarElement(prev.getLength() + length, operator));
            }
            else
            {
                cigar.add(new CigarElement(length, operator));
            }
        };

        for (int i = 0; i < Operations.length();)
        {
            char c = Operations.charAt(i);

            // Count of matching bases.
            if (isDigit(c))
            {
                int startIndex = i;
                do
                {
                    i++;
                } while (i < Operations.length() && isDigit(c = Operations.charAt(i)));
                int count = parseInt(Operations.substring(startIndex, i));
                addOperation.accept(CigarOperator.M, count);
            }

            if (i >= Operations.length())
            {
                break;
            }

            i++;
            char next = Operations.charAt(i);

            if (isLetter(c) && isLetter(next))
            {
                // Mismatch.
                addOperation.accept(CigarOperator.X, 1);
            }
            else if (isLetter(c) && next == '-')
            {
                // Insertion in query.
                addOperation.accept(CigarOperator.I, 1);
            }
            else if (c == '-' && isLetter(next))
            {
                // Deletion in query.
                addOperation.accept(CigarOperator.D, 1);
            }
            else
            {
                throw new IllegalArgumentException(format("Invalid operation: %s at %c%c", Operations, c, next));
            }
            i++;
        }
        return cigar;
    }

    public static boolean isPrimaryBlastnMatch(final BlastnMatch match)
    {
        return match.getSubjectTitle().contains("Primary Assembly") ||
                match.getSubjectTitle().contains("unlocalized genomic scaffold") ||
                match.getSubjectTitle().contains("unplaced genomic scaffold");
    }

    public static double calcSumBitScore(final Collection<BlastnMatch> matches, final int minAlignmentLength)
    {
        if(matches.isEmpty())
            return 0;

        double sumBitScore = 0;

        // process all the matches and sum up the bit score, but remove the one with best match
        Optional<BlastnMatch> bestMatch = matches.stream()
                .filter(x -> isPrimaryBlastnMatch(x))
                .max(Comparator.comparing(BlastnMatch::getBitScore));

        if(bestMatch.isPresent())
        {
            sumBitScore -= bestMatch.get().getBitScore();
        }

        for(BlastnMatch m : matches)
        {
            if(m.getAlignmentLength() >= minAlignmentLength && isPrimaryBlastnMatch(m))
            {
                sumBitScore += m.getBitScore();
            }
        }

        return sumBitScore;
    }
}
