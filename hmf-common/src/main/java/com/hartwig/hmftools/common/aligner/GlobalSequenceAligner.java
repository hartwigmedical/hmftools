package com.hartwig.hmftools.common.aligner;

import static com.hartwig.hmftools.common.aligner.WorkMatrix.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

//
// Global aligner.
// Uses the Needleman-Wunsch Algorithm to align whole sequences.
//
// Usage:
//    var aligner = new GlobalSequenceAligner(1, -3, -5, -1);
//    GlobalSequenceAligner.Alignment alignment = aligner.alignSequence(seq, refSeq);
//
// The alignment object contains a list of AlignmentOperators, which are MATCH, MISMATCH, INSERTION or DELETION.
//
// We can also use alignment.log() function to log the alignment in the format:
//
//  alignment:
//  MSMMMMMMMSMMMMMSMMMMMM-MMMMMM+MMSMMMMSMMMSMSMSMSM
//  TGCAGCTTAACTGAGAGCCGCT CTCTTCTCTTATGTTCCTTGTGTGAC
//  | ||||||| ||||| |||||| |||||| || |||| ||| | | | |
//  TCCAGCTTACCTGAGCGCCGCTCCTCTTC-CTCATGTCCCTCGGGGGCC
//
// For the scoring function, the match/mismatch score 1/-4 optimizes the scoring for 100% identical sequences and 1/-1 for 75% identical
// sequences. The default for NCBI Blastn is 2/-3, which is optimal for 89% identical sequences. BWA uses 1/-4.
// There is also gap opening and gap extension. BWA uses gap opening of -6 and gap extension of -1.
//
// Can use the following to decide good match/mismatch score for the % similarity:
// https://bioinformaticshome.com/online_software/evaluateDNAscoring/evaluateDNAscoring.html
//
public class GlobalSequenceAligner extends AlignerTraits
{
    public static class Alignment
    {
        private final String mFirstSequence;
        private final String mSecondSequence;

        private final List<AlignmentOperator> mAlignmentOperators;
        private final int mScore;

        public Alignment(String firstSequence, String secondSequence, List<AlignmentOperator> alignOps, int score)
        {
            mFirstSequence = firstSequence;
            mSecondSequence = secondSequence;
            mAlignmentOperators = alignOps;
            mScore = score;
        }

        public String getFirstSequence() { return mFirstSequence; }
        public String getSecondSequence() { return mSecondSequence; }

        public List<AlignmentOperator> getOperators() { return mAlignmentOperators; }
        public int getScore() { return mScore; }

        public String getOperatorsString() { return AlignmentOperator.toString(mAlignmentOperators); }

        public void log(@NotNull Logger logger, @NotNull Level logLevel)
        {
            AlignmentOperator.logAlignment(logger, logLevel, mFirstSequence, mSecondSequence, mAlignmentOperators);
        }
    }

    // default we use BWA settings
    // match/mismatch/gap open/gap ext of 1/-4/-6/-1
    // this score is useful for close to 100% alignment, such as when
    // aligning reads to the reference genome
    public GlobalSequenceAligner()
    {
        this(1, -4, -6, -1);
    }

    public GlobalSequenceAligner(final int matchScore, final int mismatchScore, final int gapOpeningScore, final int gapExtensionScore)
    {
        super(matchScore, mismatchScore, gapOpeningScore, gapExtensionScore);
    }

    @NotNull
    public Alignment alignSequence(@NotNull String seq, @NotNull String refSeq)
    {
        return alignSequenceImpl(seq, refSeq);
    }

    @NotNull
    private Alignment alignSequenceImpl(@NotNull String seq, @NotNull String refSeq)
    {
        int nRows = seq.length() + 1;
        int nCols = refSeq.length() + 1;
        var matrix = new WorkMatrix(nRows, nCols);

        // initialise matrix
        // first step each we initialise the first row and first column of the matrix
        matrix.setEntry(0, 0, 0, TRACEBACK_END);

        // first row, this corresponds to deletes at the start
        for (int y = 1; y < nCols; ++y)
        {
            matrix.setEntry(0, y, matrix.getScore(0, y - 1) + mGapOpeningScore, TRACEBACK_LEFT);
        }

        // first column, this corresponds to inserts at the start
        for (int x = 1; x < nRows; ++x)
        {
            matrix.setEntry(x, 0, matrix.getScore(x - 1, 0) + mGapOpeningScore, TRACEBACK_UP);
        }

        // now we can fill up the rest of the matrix as required
        for (int x = 1; x < nRows; ++x)
        {
            for (int y = 1; y < nCols; ++y)
            {
                int i = x - 1;
                int j = y - 1;

                // now we need to work out which path to take
                int diagScore = matrix.getScore(i, j);

                if (seq.charAt(i) == refSeq.charAt(j))
                {
                    diagScore += mMatchScore;
                }
                else
                {
                    diagScore += mMismatchScore;
                }

                // going left means this base is a delete
                // we must decide if this is gap extension or gap opening
                i = x;
                j = y - 1;
                int leftScore = matrix.getScore(i, j);

                if (matrix.getTraceback(i, j) == TRACEBACK_LEFT)
                {
                    // extend previous gap
                    leftScore += mGapExtensionScore;
                }
                else
                {
                    // open new gap
                    leftScore += mGapOpeningScore;
                }

                // going up means this base is an insert
                i = x - 1;
                j = y;
                int upScore = matrix.getScore(i, j);

                if (matrix.getTraceback(i, j) == TRACEBACK_UP)
                {
                    // extend previous gap
                    upScore += mGapExtensionScore;
                }
                else
                {
                    // open new gap
                    upScore += mGapOpeningScore;
                }

                //assert(!Double.isNaN(diagScore));
                //assert(!Double.isNaN(leftScore));
                //assert(!Double.isNaN(upScore));

                // for the purpose we want to
                if (diagScore >= leftScore)
                {
                    if (diagScore >= upScore)
                    {
                        matrix.setEntry(x, y, diagScore, TRACEBACK_DIAG);
                    }
                    else
                    {
                        matrix.setEntry(x, y, upScore, TRACEBACK_UP);
                    }
                }
                else if (leftScore >= upScore)
                {
                    matrix.setEntry(x, y, leftScore, TRACEBACK_LEFT);
                }
                else
                {
                    matrix.setEntry(x, y, upScore, TRACEBACK_UP);
                }

                //assert(!Double.isNaN(e.Score));
            }
        }

        var alignOps = new ArrayList<AlignmentOperator>();

        // now we apply trace back, we start from the last cell and go backwards
        for (int x = nRows - 1, y = nCols - 1; x > 0 || y > 0;)
        {
            switch (matrix.getTraceback(x, y))
            {
                case TRACEBACK_DIAG:
                    // match or substitution
                    alignOps.add(seq.charAt(x - 1) == refSeq.charAt(y - 1) ? AlignmentOperator.MATCH : AlignmentOperator.MISMATCH);
                    x--;
                    y--;
                    break;
                case TRACEBACK_LEFT:
                    // deletion
                    y--;
                    alignOps.add(AlignmentOperator.DELETION);
                    break;
                case TRACEBACK_UP:
                    // insertion
                    x--;
                    alignOps.add(AlignmentOperator.INSERTION);
                    break;
                default:
                    break;
            }
        }

        Collections.reverse(alignOps);

        int score = matrix.getScore(nRows - 1, nCols - 1);
        var alignment = new Alignment(seq, refSeq, alignOps, score);

        // log matrix
        if (mLogWorkMatrix && LOGGER.isTraceEnabled())
        {
            matrix.log(LOGGER, Level.TRACE, seq, refSeq);
        }

        if (LOGGER.isTraceEnabled())
        {
            alignment.log(LOGGER, Level.TRACE);
        }

        return alignment;
    }

    private static final Logger LOGGER = LogManager.getLogger(GlobalSequenceAligner.class);
}
