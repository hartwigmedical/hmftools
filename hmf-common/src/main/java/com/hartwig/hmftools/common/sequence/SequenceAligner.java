package com.hartwig.hmftools.common.sequence;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

/**
 * use the Needleman-Wunsch Algorithm to align sequences
 *
 * Usage: List<SequenceAligner.AlignOp> alignOps = SequenceAligner.alignSequence(seq, refSeq);
 * This gives a list of AlignOps, which are MATCH, SUBSTITUTION, INSERTION or DELETION
 * Using the align ops we can work out the alignment results.
 *
 * We can also use logAlignment function to log the alignment in the format
 *
 *
 */
public class SequenceAligner
{
    private static final Logger LOGGER = LogManager.getLogger(SequenceAligner.class);

    public enum AlignOp
    {
        MATCH('M'),
        SUBSTITUTION('S'),
        INSERTION('+'),
        DELETION('-');

        public final char code;

        AlignOp(char code) {
            this.code = code;
        }

        // convert the list of align ops to string like
        // MMMMSIDDMM
        @NotNull
        public static String toString(@NotNull List<AlignOp> alignOps)
        {
            StringBuilder b = new StringBuilder();
            alignOps.stream().map(x -> x.code).forEach(b::append);
            return b.toString();
        }
    }

    public interface IScorer
    {
        // only allow using int, this is for speed optimisation
        int scoreInsert(int seqIndex, int refSeqIndex);
        int scoreDelete(int seqIndex, int refSeqIndex);
        int scorePair(int seqIndex, int refSeqIndex, char base, char refBase);
    }

    public static class DefaultScorer implements IScorer
    {
        @Override
        public int scoreInsert(int seqIndex, int refSeqIndex)
        {
            return -1;
        }

        @Override
        public int scoreDelete(int seqIndex, int refSeqIndex)
        {
            return -1;
        }

        @Override
        public int scorePair(int seqIndex, int refSeqIndex, char base, char refBase)
        {
            if (base == refBase)
            {
                // matched pair
                return 1;
            }
            return -1;
        }
    }

    private static final int TRACEBACK_END = 0;
    private static final int TRACEBACK_DIAG = 1;
    private static final int TRACEBACK_LEFT = 2;
    private static final int TRACEBACK_UP = 3;

    private static char tracebackLabel(int move)
    {
        switch (move)
        {
            case TRACEBACK_END:
                return ' ';
            case TRACEBACK_DIAG:
                return '↖';
            case TRACEBACK_LEFT:
                return '←';
            case TRACEBACK_UP:
                return '↑';
            default:
                throw new RuntimeException("invalid traceback");
        }
    }

    // the internal work matrix used in the Needleman-Wunsch algorithm
    // NOTE: here we deliberately use a int array instead of an object
    // array to avoid slow memory access of objects allocated on the heap
    // this speeds up the algorithm by more than 100%
    private static final class Matrix
    {
        private int mNumRows;
        private int mNumCols;

        // we use a int matrix with bit shift to make it fast
        // the first 30 bits is the score, the last 2 bits is the traceback move
        private int[] mEntries;

        public Matrix(int numRows, int numCols)
        {
            mNumRows = numRows;
            mNumCols = numCols;
            mEntries = new int[numRows * numCols];
        }

        public int numRows() { return mNumRows; }
        public int numCols() { return mNumCols; }

        private int entryIndex(int row, int col)
        {
            return (row * mNumCols + col);
        }

        public int getScore(int row, int col)
        {
            return mEntries[entryIndex(row, col)] >> 2;
        }
        public int getTraceback(int row, int col)
        {
            return mEntries[entryIndex(row, col)] & 0x3;
        }
        public void setEntry(int row, int col, int score, int tracebackMove)
        {
            assert(tracebackMove == (tracebackMove & 0x3));
            assert(score == ((score << 2) >> 2));
            mEntries[entryIndex(row, col)] = (score << 2) | (tracebackMove & 0x3);
        }
    }

    @NotNull
    public static List<AlignOp> alignSequence(@NotNull String seq, @NotNull String refSeq)
    {
        return alignSequenceImpl(seq, refSeq, false, new DefaultScorer());
    }

    @NotNull
    public static List<AlignOp> alignSequence(@NotNull String seq, @NotNull String refSeq, @NotNull IScorer scorer)
    {
        return alignSequenceImpl(seq, refSeq, false, scorer);
    }

    @NotNull
    public static List<AlignOp> alignSubsequence(@NotNull String seq, @NotNull String refSeq)
    {
        return alignSequenceImpl(seq, refSeq, true, new DefaultScorer());
    }

    @NotNull
    public static List<AlignOp> alignSubsequence(@NotNull String seq, @NotNull String refSeq, @NotNull IScorer scorer)
    {
        return alignSequenceImpl(seq, refSeq, true, scorer);
    }

    @NotNull
    private static List<AlignOp> alignSequenceImpl(@NotNull String seq, @NotNull String refSeq, boolean isSubsequence, @NotNull IScorer scorer)
    {
        int nRows = seq.length() + 1;
        int nCols = refSeq.length() + 1;
        Matrix matrix = new Matrix(nRows, nCols);

        // initialise matrix
        // first step each we initialise the first row and first column of the matrix
        matrix.setEntry(0, 0, 0, TRACEBACK_END);

        // first row, this corresponds to deletes at the start
        for (int y = 1; y < nCols; ++y)
        {
            // Matrix.Entry e = matrix.at(0, y);

            if (isSubsequence)
            {
                // NOTE: special case:
                // we want to match it such that the ref seq end can extend beyond the end at no cost
                // i.e.
                // ---TTAGACGTC
                //    |||||||||
                // CGTTTAGACGTC
                // we want the 3 deletions at the start to be free
                // in the standard alogrithm, we would use
                matrix.setEntry(0, y, matrix.getScore(0, y - 1), TRACEBACK_LEFT);
            }
            else
            {
                matrix.setEntry(0, y, matrix.getScore(0, y - 1) + scorer.scoreDelete(0, y - 1), TRACEBACK_LEFT);
            }
        }

        // first column, this corresponds to inserts at the start
        for (int x = 1; x < nRows; ++x)
        {
            matrix.setEntry(x, 0, matrix.getScore(x - 1, 0) + scorer.scoreInsert(x - 1, 0), TRACEBACK_UP);
        }

        // now we can fill up the rest of the matrix as required
        for (int x = 1; x < nRows; ++x)
        {
            for (int y = 1; y < nCols; ++y)
            {
                int i = x - 1;
                int j = y - 1;

                // now we need to work out which path to take
                int diagScore = (matrix.getScore(i, j) + scorer.scorePair(i, j, seq.charAt(i), refSeq.charAt(j)));

                i = x;
                j = y - 1;
                int leftScore;
                if (isSubsequence && i == (nRows - 1))
                {
                    // NOTE: special case:
                    // we want to match it such that the ref seq end can extend beyond the end at no cost
                    // i.e.
                    // TTAGACGTC----
                    // |||||||||
                    // TTAGACGTCCGTC
                    // we want the 4 deletions at the end to be free
                    // reason is that we are trying to align to a telomere template, and we would provide
                    // a template that is likely to be extended in both ends
                    leftScore = matrix.getScore(i, j);
                }
                else
                {
                    leftScore = matrix.getScore(i, j) + scorer.scoreDelete(i, j);
                }

                i = x - 1;
                j = y;
                int upScore = matrix.getScore(i, j) + scorer.scoreInsert(i, j);

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

        List<AlignOp> alignOps = new ArrayList<>();

        // now we apply trace back, we start from the last cell and go backwards
        for (int x = nRows - 1, y = nCols - 1; x > 0 || y > 0;)
        {
            switch (matrix.getTraceback(x, y))
            {
                case TRACEBACK_DIAG:
                    // match or substitution
                    if (seq.charAt(x - 1) == refSeq.charAt(y - 1))
                    {
                        alignOps.add(AlignOp.MATCH);
                    } else
                    {
                        alignOps.add(AlignOp.SUBSTITUTION);
                    }
                    x--;
                    y--;
                    break;
                case TRACEBACK_LEFT:
                    // deletion
                    y--;
                    alignOps.add(AlignOp.DELETION);
                    break;
                case TRACEBACK_UP:
                    // insertion
                    x--;
                    alignOps.add(AlignOp.INSERTION);
                    break;
                case TRACEBACK_END:
                    break;
            }
        }

        Collections.reverse(alignOps);

        // log matrix
        /* if (LOGGER.isTraceEnabled())
        {
            logWorkMatrix(Level.TRACE, seq, refSeq, matrix);
        }*/

        if (LOGGER.isTraceEnabled())
        {
            logAlignment(Level.TRACE, seq, refSeq, alignOps);
        }

        return alignOps;
    }

    public static void logWorkMatrix(@NotNull Level logLevel, @NotNull String seq, @NotNull String refSeq, @NotNull Matrix matrix)
    {
        LOGGER.log(logLevel, "alignment matrix:");

        StringBuilder lineBuilder = new StringBuilder();
        lineBuilder.append("       ");

        for (int y = 1; y < matrix.numCols(); ++y)
        {
            lineBuilder.append(String.format("  %c  ", refSeq.charAt(y - 1)));
        }

        LOGGER.log(logLevel, "{}", lineBuilder);

        for (int x = 0; x < matrix.numRows(); ++x)
        {
            lineBuilder = new StringBuilder();
            if (x >= 1)
                lineBuilder.append(String.format("%c ", seq.charAt(x - 1)));
            else
                lineBuilder.append("  ");

            for (int y = 0; y < matrix.numCols(); ++y)
            {
                lineBuilder.append(String.format("%3d%c ", matrix.getScore(x, y), tracebackLabel(matrix.getTraceback(x, y))));
            }

            LOGGER.log(logLevel, "{}", lineBuilder);
        }
    }

    public static void logAlignment(@NotNull Level logLevel, @NotNull String seq, @NotNull String template, @NotNull List<AlignOp> alignOps)
    {
        if (!LOGGER.isEnabled(logLevel))
            return;

        StringBuilder seqAlignBuilder = new StringBuilder();
        StringBuilder alignOpBuilder = new StringBuilder();
        StringBuilder refSeqAlignBuilder = new StringBuilder();

        int i = 0;
        int j = 0;

        for (AlignOp op : alignOps)
        {
            switch (op)
            {
                case MATCH:
                    seqAlignBuilder.append(seq.charAt(i++));
                    refSeqAlignBuilder.append(template.charAt(j++));
                    alignOpBuilder.append('|');
                    break;
                case SUBSTITUTION:
                    seqAlignBuilder.append(seq.charAt(i++));
                    refSeqAlignBuilder.append(template.charAt(j++));
                    alignOpBuilder.append(' ');
                    break;
                case INSERTION:
                    seqAlignBuilder.append(seq.charAt(i++));
                    refSeqAlignBuilder.append('-');
                    alignOpBuilder.append(' ');
                    break;
                case DELETION:
                    seqAlignBuilder.append('-');
                    refSeqAlignBuilder.append(template.charAt(j++));
                    alignOpBuilder.append(' ');
                    break;
            }
        }

        LOGGER.log(logLevel, "alignment:");
        LOGGER.log(logLevel, "{}", AlignOp.toString(alignOps));
        LOGGER.log(logLevel, "{}", seqAlignBuilder);
        LOGGER.log(logLevel, "{}", alignOpBuilder);
        LOGGER.log(logLevel, "{}", refSeqAlignBuilder);
    }
}
