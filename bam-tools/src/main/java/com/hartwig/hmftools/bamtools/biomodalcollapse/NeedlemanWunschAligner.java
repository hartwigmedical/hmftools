package com.hartwig.hmftools.bamtools.biomodalcollapse;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.M;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.PREFIX_X;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.PREFIX_Y;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.SUFFIX_X;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.SUFFIX_Y;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.X;
import static com.hartwig.hmftools.bamtools.biomodalcollapse.NeedlemanWunschAligner.TraceBackState.Y;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.OptionalInt;
import java.util.function.ToIntBiFunction;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class NeedlemanWunschAligner<T>
{
    private long[] MScore;
    private long[] XScore;
    private long[] YScore;
    private long[] PreXScore;
    private long[] PreYScore;
    private long[] SuffXScore;
    private long[] SuffYScore;

    private TraceBackState[] MTrace;
    private TraceBackState[] XTrace;
    private TraceBackState[] YTrace;
    private TraceBackState[] PreXTrace;
    private TraceBackState[] PreYTrace;
    private TraceBackState[] SuffXTrace;
    private TraceBackState[] SuffYTrace;

    public NeedlemanWunschAligner()
    {
        MScore = null;
        XScore = null;
        YScore = null;
        PreXScore = null;
        PreYScore = null;
        SuffXScore = null;
        SuffYScore = null;

        MTrace = null;
        XTrace = null;
        YTrace = null;
        PreXTrace = null;
        PreYTrace = null;
        SuffXTrace = null;
        SuffYTrace = null;
    }

    public enum TraceBackState
    {
        M,
        X,
        Y,
        PREFIX_X,
        PREFIX_Y,
        SUFFIX_X,
        SUFFIX_Y
    }

    private void initializeMatrices(int requiredSize)
    {
        if(MScore != null && MScore.length >= requiredSize)
        {
            return;
        }

        int size = MScore == null ? 1 : MScore.length;
        while(size < requiredSize)
        {
            size *= 2;
        }

        MScore = new long[size];
        XScore = new long[size];
        YScore = new long[size];
        PreXScore = new long[size];
        PreYScore = new long[size];
        SuffXScore = new long[size];
        SuffYScore = new long[size];

        MTrace = new TraceBackState[size];
        XTrace = new TraceBackState[size];
        YTrace = new TraceBackState[size];
        PreXTrace = new TraceBackState[size];
        PreYTrace = new TraceBackState[size];
        SuffXTrace = new TraceBackState[size];
        SuffYTrace = new TraceBackState[size];
    }

    @Nullable
    public List<Pair<T, T>> align(final List<T> seq1, final List<T> seq2, final ToIntBiFunction<T, T> scoringFn, int openGapPenalty,
            int extendGapPenalty, boolean noPenaltySeq1Prefix, boolean noPenaltySeq2Prefix, boolean noPenaltySeq1Suffix,
            boolean noPenaltySeq2Suffix, final OptionalInt maxIndelBandWidth)
    {
        initializeMatrices((seq1.size() + 1) * (seq2.size() + 1));

        int minDiagIndex = 0;
        int maxDiagIndex = 0;
        if(maxIndelBandWidth.isPresent())
        {
            minDiagIndex = min(0, seq1.size() - seq2.size()) - maxIndelBandWidth.getAsInt();
            maxDiagIndex = max(0, seq1.size() - seq2.size()) + maxIndelBandWidth.getAsInt();
        }

        int minRow = maxIndelBandWidth.isPresent() ? max(1, minDiagIndex) : 1;
        int maxRow = maxIndelBandWidth.isPresent() ? min(seq1.size(), maxDiagIndex) : seq1.size();
        for(int r = minRow; r <= maxRow; r++)
        {
            int index = (seq2.size() + 1) * r;

            MScore[index] = Integer.MIN_VALUE;
            XScore[index] = -openGapPenalty - (long) (r - 1) * extendGapPenalty;
            YScore[index] = Integer.MIN_VALUE;

            MTrace[index] = null;
            XTrace[index] = (r == 1) ? M : X;
            YTrace[index] = null;

            if(noPenaltySeq1Prefix)
            {
                PreXScore[index] = 0;
                PreXTrace[index] = (r == 1) ? M : PREFIX_X;
            }

            if(noPenaltySeq2Prefix)
            {
                PreYScore[index] = Integer.MIN_VALUE;
                PreYTrace[index] = null;
            }

            if(noPenaltySeq1Suffix)
            {
                SuffXScore[index] = 0;
                SuffXTrace[index] = (r == 1) ? M : SUFFIX_X;
            }

            if(noPenaltySeq2Suffix)
            {
                SuffYScore[index] = Integer.MIN_VALUE;
                SuffYTrace[index] = null;
            }
        }

        int minCol = maxIndelBandWidth.isPresent() ? max(1, -maxDiagIndex) : 1;
        int maxCol = maxIndelBandWidth.isPresent() ? min(seq2.size(), -minDiagIndex) : seq2.size();
        for(int c = minCol; c <= maxCol; c++)
        {
            MScore[c] = Integer.MIN_VALUE;
            XScore[c] = Integer.MIN_VALUE;
            YScore[c] = -openGapPenalty - (long) (c - 1) * extendGapPenalty;

            MTrace[c] = null;
            XTrace[c] = null;
            YTrace[c] = (c == 1) ? M : Y;

            if(noPenaltySeq1Prefix)
            {
                PreXScore[c] = Integer.MIN_VALUE;
                PreXTrace[c] = null;
            }

            if(noPenaltySeq2Prefix)
            {
                PreYScore[c] = 0;
                PreYTrace[c] = (c == 1) ? M : PREFIX_Y;
            }

            if(noPenaltySeq1Suffix)
            {
                SuffXScore[c] = Integer.MIN_VALUE;
                SuffXTrace[c] = null;
            }

            if(noPenaltySeq2Suffix)
            {
                SuffYScore[c] = 0;
                SuffYTrace[c] = (c == 1) ? M : SUFFIX_Y;
            }
        }

        MScore[0] = 0;
        XScore[0] = Integer.MIN_VALUE;
        YScore[0] = Integer.MIN_VALUE;

        MTrace[0] = null;
        XTrace[0] = null;
        YTrace[0] = null;

        if(noPenaltySeq1Prefix)
        {
            PreXScore[0] = Integer.MIN_VALUE;
            PreXTrace[0] = null;
        }

        if(noPenaltySeq2Prefix)
        {
            PreYScore[0] = Integer.MIN_VALUE;
            PreYTrace[0] = null;
        }

        if(noPenaltySeq1Suffix)
        {
            SuffXScore[0] = Integer.MIN_VALUE;
            SuffXTrace[0] = null;
        }

        if(noPenaltySeq2Suffix)
        {
            SuffYScore[0] = Integer.MIN_VALUE;
            SuffYTrace[0] = null;
        }

        // fill out
        List<Pair<TraceBackState, Long>> edgeScores = Lists.newArrayList();
        for(int r = 1; r <= seq1.size(); r++)
        {
            T base1 = seq1.get(r - 1);

            minCol = maxIndelBandWidth.isPresent() ? max(1, -maxDiagIndex + r) : 1;
            maxCol = maxIndelBandWidth.isPresent() ? min(seq2.size(), -minDiagIndex + r) : seq2.size();
            for(int c = minCol; c <= maxCol; c++)
            {
                int diag = r - c;

                T base2 = seq2.get(c - 1);
                int index = (seq2.size() + 1) * r + c;
                int upIndex = diag == minDiagIndex && maxIndelBandWidth.isPresent() ? -1 : (seq2.size() + 1) * (r - 1) + c;
                int leftIndex = diag == maxDiagIndex && maxIndelBandWidth.isPresent() ? -1 : (seq2.size() + 1) * r + c - 1;
                int diagIndex = (seq2.size() + 1) * (r - 1) + c - 1;

                // M
                int matchScore = scoringFn.applyAsInt(base1, base2);
                edgeScores.clear();
                edgeScores.add(Pair.of(M, MScore[diagIndex] + matchScore));
                edgeScores.add(Pair.of(X, XScore[diagIndex] + matchScore));
                edgeScores.add(Pair.of(Y, YScore[diagIndex] + matchScore));
                edgeScores.add(Pair.of(PREFIX_X, noPenaltySeq1Prefix ? PreXScore[diagIndex] + matchScore : Integer.MIN_VALUE));
                edgeScores.add(Pair.of(PREFIX_Y, noPenaltySeq2Prefix ? PreYScore[diagIndex] + matchScore : Integer.MIN_VALUE));
                Pair<TraceBackState, Long> bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                MScore[index] = bestTransition.getRight();
                MTrace[index] = bestTransition.getLeft();

                // X
                edgeScores.clear();
                edgeScores.add(Pair.of(M, upIndex >= 0 ? MScore[upIndex] - openGapPenalty : Integer.MIN_VALUE));
                edgeScores.add(Pair.of(X, upIndex >= 0 ? XScore[upIndex] - extendGapPenalty : Integer.MIN_VALUE));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                XScore[index] = bestTransition.getRight();
                XTrace[index] = bestTransition.getLeft();

                // Y
                edgeScores.clear();
                edgeScores.add(Pair.of(M, leftIndex >= 0 ? MScore[leftIndex] - openGapPenalty : Integer.MIN_VALUE));
                edgeScores.add(Pair.of(Y, leftIndex >= 0 ? YScore[leftIndex] - extendGapPenalty : Integer.MIN_VALUE));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                YScore[index] = bestTransition.getRight();
                YTrace[index] = bestTransition.getLeft();

                if(noPenaltySeq1Prefix)
                {
                    // preX
                    PreXScore[index] = upIndex >= 0 ? PreXScore[upIndex] : Integer.MIN_VALUE;
                    PreXTrace[index] = PREFIX_X;
                }

                if(noPenaltySeq2Prefix)
                {
                    // preY
                    PreYScore[index] = leftIndex >= 0 ? PreYScore[leftIndex] : Integer.MIN_VALUE;
                    PreYTrace[index] = PREFIX_Y;
                }

                if(noPenaltySeq1Suffix)
                {
                    // suffX
                    edgeScores.clear();
                    edgeScores.add(Pair.of(M, upIndex >= 0 ? MScore[upIndex] : Integer.MIN_VALUE));
                    edgeScores.add(Pair.of(SUFFIX_X, upIndex >= 0 ? SuffXScore[upIndex] : Integer.MIN_VALUE));
                    edgeScores.add(Pair.of(PREFIX_Y, noPenaltySeq2Prefix && upIndex >= 0 ? PreYScore[upIndex] : Integer.MIN_VALUE));
                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                    SuffXScore[index] = bestTransition.getRight();
                    SuffXTrace[index] = bestTransition.getLeft();
                }

                if(noPenaltySeq2Suffix)
                {
                    // suffY
                    edgeScores.clear();
                    edgeScores.add(Pair.of(M, leftIndex >= 0 ? MScore[leftIndex] : Integer.MIN_VALUE));
                    edgeScores.add(Pair.of(SUFFIX_Y, leftIndex >= 0 ? SuffYScore[leftIndex] : Integer.MIN_VALUE));
                    edgeScores.add(Pair.of(PREFIX_X, noPenaltySeq1Prefix && leftIndex >= 0 ? PreXScore[leftIndex] : Integer.MIN_VALUE));
                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                    SuffYScore[index] = bestTransition.getRight();
                    SuffYTrace[index] = bestTransition.getLeft();
                }
            }
        }

        // trace back
        List<Pair<T, T>> alignment = Lists.newArrayList();
        int r = seq1.size();
        int c = seq2.size();
        int index = (seq2.size() + 1) * r + c;

        List<Pair<TraceBackState, Long>> finalScores = Lists.newArrayList();
        finalScores.add(Pair.of(M, MScore[index]));
        finalScores.add(Pair.of(X, XScore[index]));
        finalScores.add(Pair.of(Y, YScore[index]));
        finalScores.add(Pair.of(PREFIX_X, noPenaltySeq1Prefix ? PreXScore[index] : Integer.MIN_VALUE));
        finalScores.add(Pair.of(PREFIX_Y, noPenaltySeq2Prefix ? PreYScore[index] : Integer.MIN_VALUE));
        finalScores.add(Pair.of(SUFFIX_X, noPenaltySeq1Suffix ? SuffXScore[index] : Integer.MIN_VALUE));
        finalScores.add(Pair.of(SUFFIX_Y, noPenaltySeq2Suffix ? SuffYScore[index] : Integer.MIN_VALUE));

        Pair<TraceBackState, Long> bestFinal = Collections.max(finalScores, Comparator.comparingLong(x -> x.getRight()));
        TraceBackState state = bestFinal.getLeft();
        while(r != 0 || c != 0)
        {
            int diag = r - c;
            if(maxIndelBandWidth.isPresent() && (diag < minDiagIndex || diag > maxDiagIndex))
            {
                throw new RuntimeException("Traceback exited prescribed band width");
            }

            if(maxIndelBandWidth.isPresent() && (diag == minDiagIndex || diag == maxDiagIndex))
            {
                return null;
            }

            index = (seq2.size() + 1) * r + c;
            if(state == M)
            {
                alignment.add(Pair.of(seq1.get(r - 1), seq2.get(c - 1)));
                state = MTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
                c--;
            }
            else if(state == X)
            {
                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = XTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == Y)
            {
                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = YTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                c--;
            }
            else if(state == PREFIX_X)
            {
                if(!noPenaltySeq1Prefix)
                {
                    throw new RuntimeException("PREFIX_X is an invalid traceback state");
                }

                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = PreXTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == PREFIX_Y)
            {
                if(!noPenaltySeq2Prefix)
                {
                    throw new RuntimeException("PREFIX_Y is an invalid traceback state");
                }

                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = PreYTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                c--;
            }
            else if(state == SUFFIX_X)
            {
                if(!noPenaltySeq1Suffix)
                {
                    throw new RuntimeException("SUFFIX_X is an invalid traceback state");
                }

                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = SuffXTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == SUFFIX_Y)
            {
                if(!noPenaltySeq2Suffix)
                {
                    throw new RuntimeException("SUFFIX_Y is an invalid traceback state");
                }

                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = SuffYTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                c--;
            }
            else
            {
                throw new RuntimeException("Unreachable");
            }
        }

        Collections.reverse(alignment);
        return alignment;
    }

    public List<Pair<T, T>> approxAlign(final List<T> seq1, final List<T> seq2, final ToIntBiFunction<T, T> scoringFn, int openGapPenalty,
            int extendGapPenalty, boolean noPenaltySeq1Prefix, boolean noPenaltySeq2Prefix, boolean noPenaltySeq1Suffix,
            boolean noPenaltySeq2Suffix)
    {
        int maxIndelBandWidth = 1;
        List<Pair<T, T>> alignment = null;
        while(alignment == null)
        {
            alignment = align(seq1, seq2, scoringFn, openGapPenalty, extendGapPenalty, noPenaltySeq1Prefix, noPenaltySeq2Prefix,
                    noPenaltySeq1Suffix, noPenaltySeq2Suffix, OptionalInt.of(maxIndelBandWidth));

            maxIndelBandWidth *= 2;
        }

        return alignment;
    }
}
