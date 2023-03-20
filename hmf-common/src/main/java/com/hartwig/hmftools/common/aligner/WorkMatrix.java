package com.hartwig.hmftools.common.aligner;

import org.apache.logging.log4j.Level;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

// the internal work matrix used in the Needleman-Wunsch and Smith-Waterman algorithms
// NOTE: here we deliberately use a int array instead of an object
// array to avoid slow memory access of objects allocated on the heap
// this speeds up the algorithm by more than 100%
class WorkMatrix
{
    static final int TRACEBACK_END = 0;
    static final int TRACEBACK_DIAG = 1;
    static final int TRACEBACK_LEFT = 2;
    static final int TRACEBACK_UP = 3;

    static char tracebackLabel(int move)
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

    private final int mNumRows;
    private final int mNumCols;

    // we use a int matrix with bit shift to make it fast
    // the first 30 bits is the score, the last 2 bits is the traceback move
    private final int[] mEntries;

    public WorkMatrix(int numRows, int numCols)
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

    public void log(@NotNull Logger logger, @NotNull Level logLevel, @NotNull String seq, @NotNull String refSeq)
    {
        if (!logger.isEnabled(logLevel))
            return;

        logger.log(logLevel, "alignment matrix:");

        StringBuilder lineBuilder = new StringBuilder();
        lineBuilder.append("       ");

        for (int y = 1; y < numCols(); ++y)
        {
            lineBuilder.append(String.format("  %c  ", refSeq.charAt(y - 1)));
        }

        logger.log(logLevel, "{}", lineBuilder);

        for (int x = 0; x < numRows(); ++x)
        {
            lineBuilder = new StringBuilder();
            if (x >= 1)
                lineBuilder.append(String.format("%c ", seq.charAt(x - 1)));
            else
                lineBuilder.append("  ");

            for (int y = 0; y < numCols(); ++y)
            {
                lineBuilder.append(String.format("%3d%c ", getScore(x, y), tracebackLabel(getTraceback(x, y))));
            }

            logger.log(logLevel, "{}", lineBuilder);
        }
    }
}
