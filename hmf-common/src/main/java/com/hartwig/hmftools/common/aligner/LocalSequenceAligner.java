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
// Find local alignments. Uses the Smith-Waterman Algorithm.
//
// Usage:
//    var aligner = new LocalSequenceAligner(1, -3, -5, -1);
//    LocalSequenceAligner.Alignment alignment = aligner.alignSequence(seq, refSeq);
//
// The alignment object contains a list of AlignmentOperators, which are MATCH, MISMATCH, INSERTION or DELETION. It also contains
// the start and end indices of the aligned segments.
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
// A good interactive tool to understand Smith-Waterman can be found here:
// https://cse442-17f.github.io/Prims-and-A-Star/
//
public class LocalSequenceAligner extends AlignerTraits
{
    public static class Alignment
    {
        private final String mFirstSequence;
        private final String mSecondSequence;

        private final int mFirstSeqAlignStart;
        private final int mFirstSeqAlignEnd;
        private final int mSecondSeqAlignStart;
        private final int mSecondSeqAlignEnd;

        private final List<AlignmentOperator> mAlignmentOperators;
        private final int mScore;

        public Alignment(
                String firstSequence, String secondSequence,
                int leftSeqAlignStart, int leftSeqAlignEnd,
                int rightSeqAlignStart, int rightSeqAlignEnd,
                List<AlignmentOperator> alignOps, int score)
        {
            mFirstSequence = firstSequence;
            mSecondSequence = secondSequence;
            mFirstSeqAlignStart = leftSeqAlignStart;
            mFirstSeqAlignEnd = leftSeqAlignEnd;
            mSecondSeqAlignStart = rightSeqAlignStart;
            mSecondSeqAlignEnd = rightSeqAlignEnd;
            mAlignmentOperators = alignOps;
            mScore = score;
        }

        public String getFirstSequence() { return mFirstSequence; }
        public String getSecondSequence() { return mSecondSequence; }

        public int getFirstSequenceAlignStart() { return mFirstSeqAlignStart; }
        public int getFirstSequenceAlignEnd() { return mFirstSeqAlignEnd; }
        public int getSecondSequenceAlignStart() { return mSecondSeqAlignStart; }
        public int getSecondSequenceAlignEnd() { return mSecondSeqAlignEnd; }

        public int getFirstSequenceAlignLength() { return mFirstSeqAlignEnd - mFirstSeqAlignStart; }
        public int getSecondSequenceAlignLength() { return mSecondSeqAlignEnd - mSecondSeqAlignStart; }

        public List<AlignmentOperator> getOperators() { return mAlignmentOperators; }
        public int getScore() { return mScore; }

        public String getOperatorsString() { return AlignmentOperator.toString(mAlignmentOperators); }

        public void log(@NotNull Logger logger, @NotNull Level logLevel)
        {
            AlignmentOperator.logAlignment(logger, logLevel, mFirstSequence, mSecondSequence,
                    mFirstSeqAlignStart, mSecondSeqAlignStart, mAlignmentOperators);
        }
    }

    // default we use BWA settings
    // match/mismatch/gap open/gap ext of 1/-4/-6/-1
    // this score is useful for close to 100% alignment, such as when
    // aligning reads to the reference genome
    public LocalSequenceAligner()
    {
        this(1, -4, -6, -1);
    }

    public LocalSequenceAligner(final int matchScore, final int mismatchScore, final int gapOpeningScore, final int gapExtensionScore)
    {
        super(matchScore, mismatchScore, gapOpeningScore, gapExtensionScore);
    }

    @NotNull
    public Alignment alignSequence(@NotNull String seq1, @NotNull String seq2)
    {
        return alignSequenceImpl(seq1, seq2);
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

        // first row all 0
        for (int y = 1; y < nCols; ++y)
        {
            matrix.setEntry(0, y, matrix.getScore(0, y - 1), TRACEBACK_END);
        }

        // first column all 0
        for (int x = 1; x < nRows; ++x)
        {
            matrix.setEntry(x, 0, matrix.getScore(x - 1, 0), TRACEBACK_END);
        }

        int highestScoreX = 0;
        int highestScoreY = 0;
        int highestScore = 0;

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

                // for S-W aligner, there are 4 choices
                // 1. diag match
                // 2. left delete
                // 3. up insert
                // 4. 0 means trace back end
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

                int score = matrix.getScore(x, y);

                // negative score means trace back end, and set to 0
                if (score <= 0)
                {
                    matrix.setEntry(x, y, 0, TRACEBACK_END);
                }
                else if (score > highestScore)
                {
                    highestScore = score;
                    highestScoreX = x;
                    highestScoreY = y;
                }
                //assert(!Double.isNaN(e.Score));
            }
        }

        var alignOps = new ArrayList<AlignmentOperator>();

        int seq1AlignStart = -1;
        int seq2AlignStart = -1;

        // now we apply trace back, we start from the cell with highest score
        int x, y;
        for (x = highestScoreX, y = highestScoreY; x >= 0 && y >= 0 && seq1AlignStart == -1;)
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
                case TRACEBACK_END:
                    // this is where the alignment starts
                    seq1AlignStart = x;
                    seq2AlignStart = y;
                    break;
            }
        }

        Collections.reverse(alignOps);

        @SuppressWarnings("SuspiciousNameCombination")
        var alignment = new Alignment(seq, refSeq, seq1AlignStart, highestScoreX, seq2AlignStart, highestScoreY, alignOps, highestScore);

        // log matrix
        if (mLogWorkMatrix && LOGGER.isTraceEnabled())
        {
            matrix.log(LOGGER, Level.TRACE, seq, refSeq);
        }

        if (LOGGER.isTraceEnabled())
        {
            alignment.log(LOGGER, Level.TRACE);
            LOGGER.trace("score: {}", alignment.getScore());
        }

        return alignment;
    }

    private static final Logger LOGGER = LogManager.getLogger(LocalSequenceAligner.class);
}
