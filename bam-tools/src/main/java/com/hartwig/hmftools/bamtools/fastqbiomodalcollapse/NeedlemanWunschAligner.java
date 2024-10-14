package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.ToIntBiFunction;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.tuple.Pair;

public class NeedlemanWunschAligner<T>
{
    private final ToIntBiFunction<T, T> mScoringFn;
    private final ToIntBiFunction<T, T> mScoringFn2;
    private final int mOpenGapPenalty;
    private final int mExtendGapPenalty;

    public NeedlemanWunschAligner(final ToIntBiFunction<T, T> scoringFn, final ToIntBiFunction<T, T> scoringFn2, int openGapPenalty, int extendGapPenalty)
    {
        mScoringFn = scoringFn;
        mScoringFn2 = scoringFn2;
        mOpenGapPenalty = openGapPenalty;
        mExtendGapPenalty = extendGapPenalty;
    }

    public enum TraceBackState
    {
        M,
        X,
        Y,
        M2,
        X2,
        Y2,
        SUFFIX_X,
        SUFFIX_Y;
    }


    public List<Pair<T, T>> align(final List<T> seq1, final List<T> seq2, boolean noPenaltySeq1Suffix, boolean noPenaltySeq2Suffix)
    {
        // initialize scoring matrices
        long[] M = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] X = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] Y = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] M2 = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] X2 = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] Y2 = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] suffX = noPenaltySeq1Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        long[] suffY = noPenaltySeq2Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        TraceBackState[] MTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] XTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] YTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] M2Trace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] X2Trace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] Y2Trace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] suffXTrace = noPenaltySeq1Suffix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        TraceBackState[] suffYTrace = noPenaltySeq2Suffix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        for(int r = 1; r <= seq1.size(); r++)
        {
            int index = (seq2.size() + 1) * r;

            M[index] = Integer.MIN_VALUE;
            X[index] = -mOpenGapPenalty - (r - 1) * mExtendGapPenalty;
            Y[index] = Integer.MIN_VALUE;
            M2[index] = Integer.MIN_VALUE;
            X2[index] = Integer.MIN_VALUE;
            Y2[index] = Integer.MIN_VALUE;

            MTrace[index] = null;
            XTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.X;
            YTrace[index] = null;
            M2Trace[index] = null;
            X2Trace[index] = null;
            Y2Trace[index] = null;

            if(noPenaltySeq1Suffix)
            {
                suffX[index] = 0;
                suffXTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.SUFFIX_X;
            }

            if(noPenaltySeq2Suffix)
            {
                suffY[index] = Integer.MIN_VALUE;
                suffYTrace[index] = null;
            }
        }

        for(int c = 1; c <= seq2.size(); c++)
        {
            M[c] = Integer.MIN_VALUE;
            X[c] = Integer.MIN_VALUE;
            Y[c] = -mOpenGapPenalty - (c - 1) * mExtendGapPenalty;
            M2[c] = Integer.MIN_VALUE;
            X2[c] = Integer.MIN_VALUE;
            Y2[c] = Integer.MIN_VALUE;

            MTrace[c] = null;
            XTrace[c] = null;
            YTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.Y;
            M2Trace[c] = null;
            X2Trace[c] = null;
            Y2Trace[c] = null;

            if(noPenaltySeq1Suffix)
            {
                suffX[c] = Integer.MIN_VALUE;
                suffXTrace[c] = null;
            }

            if(noPenaltySeq2Suffix)
            {
                suffY[c] = 0;
                suffYTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.SUFFIX_Y;
            }
        }

        M[0] = 0;
        X[0] = Integer.MIN_VALUE;
        Y[0] = Integer.MIN_VALUE;
        M2[0] = Integer.MIN_VALUE;
        X2[0] = Integer.MIN_VALUE;
        Y2[0] = Integer.MIN_VALUE;

        MTrace[0] = null;
        XTrace[0] = null;
        YTrace[0] = null;
        M2Trace[0] = null;
        X2Trace[0] = null;
        Y2Trace[0] = null;

        if(noPenaltySeq1Suffix)
        {
            suffX[0] = Integer.MIN_VALUE;
            suffXTrace[0] = null;
        }

        if(noPenaltySeq2Suffix)
        {
            suffY[0] = Integer.MIN_VALUE;
            suffYTrace[0] = null;
        }

        // fill out
        for(int r = 1; r <= seq1.size(); r++)
        {
            T base1 = seq1.get(r - 1);
            for(int c = 1; c <= seq2.size(); c++)
            {
                T base2 = seq2.get(c - 1);
                int index = (seq2.size() + 1) * r + c;
                int upIndex = (seq2.size() + 1) * (r - 1) + c;
                int leftIndex = (seq2.size() + 1) * r + c - 1;
                int diagIndex = (seq2.size() + 1) * (r - 1) + c - 1;
                int matchScore = mScoringFn.applyAsInt(base1, base2);
                int matchScore2 = mScoringFn2.applyAsInt(base1, base2);

                List<Pair<TraceBackState, Long>> edgeScores = Lists.newArrayList();

                // M
                edgeScores.add(Pair.of(TraceBackState.M, M[diagIndex] + matchScore));
                edgeScores.add(Pair.of(TraceBackState.X, X[diagIndex] + matchScore));
                edgeScores.add(Pair.of(TraceBackState.Y, Y[diagIndex] + matchScore));
                Pair<TraceBackState, Long> bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                M[index] = bestTransition.getRight();
                MTrace[index] = bestTransition.getLeft();

                // X
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[upIndex] - mOpenGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.X, X[upIndex] - mExtendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                X[index] = bestTransition.getRight();
                XTrace[index] = bestTransition.getLeft();

                // Y
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex] - mOpenGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.Y, Y[leftIndex] - mExtendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                Y[index] = bestTransition.getRight();
                YTrace[index] = bestTransition.getLeft();

                // M2
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[diagIndex] + matchScore2));
                edgeScores.add(Pair.of(TraceBackState.X, X[diagIndex] + matchScore2));
                edgeScores.add(Pair.of(TraceBackState.Y, Y[diagIndex] + matchScore2));
                edgeScores.add(Pair.of(TraceBackState.M2, M2[diagIndex] + matchScore2));
                edgeScores.add(Pair.of(TraceBackState.X2, X2[diagIndex] + matchScore2));
                edgeScores.add(Pair.of(TraceBackState.Y2, Y2[diagIndex] + matchScore2));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                M2[index] = bestTransition.getRight();
                M2Trace[index] = bestTransition.getLeft();

                // X2
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M2, M2[upIndex] - mOpenGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.X2, X2[upIndex] - mExtendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                X2[index] = bestTransition.getRight();
                X2Trace[index] = bestTransition.getLeft();

                // Y2
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M2, M2[leftIndex] - mOpenGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.Y2, Y2[leftIndex] - mExtendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                Y2[index] = bestTransition.getRight();
                Y2Trace[index] = bestTransition.getLeft();

                // suffX
                if(noPenaltySeq1Suffix)
                {
                    edgeScores.clear();
                    edgeScores.add(Pair.of(TraceBackState.M, M[upIndex]));
                    edgeScores.add(Pair.of(TraceBackState.M2, M2[upIndex]));
                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_X, suffX[upIndex]));
                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                    suffX[index] = bestTransition.getRight();
                    suffXTrace[index] = bestTransition.getLeft();
                }

                // suffY
                if(noPenaltySeq2Suffix)
                {
                    edgeScores.clear();
                    edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex]));
                    edgeScores.add(Pair.of(TraceBackState.M2, M2[leftIndex]));
                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_Y, suffY[leftIndex]));
                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                    suffY[index] = bestTransition.getRight();
                    suffYTrace[index] = bestTransition.getLeft();
                }
            }
        }

        // trace back
        List<Pair<T, T>> alignment = Lists.newArrayList();
        int r = seq1.size();
        int c = seq2.size();
        int index = (seq2.size() + 1) * r + c;

        List<Pair<TraceBackState, Long>> finalScores = Lists.newArrayList();
        finalScores.add(Pair.of(TraceBackState.M, M[index]));
        finalScores.add(Pair.of(TraceBackState.X, X[index]));
        finalScores.add(Pair.of(TraceBackState.Y, Y[index]));
        finalScores.add(Pair.of(TraceBackState.M2, M2[index]));
        finalScores.add(Pair.of(TraceBackState.X2, X2[index]));
        finalScores.add(Pair.of(TraceBackState.Y2, Y2[index]));
        finalScores.add(Pair.of(TraceBackState.SUFFIX_X, noPenaltySeq1Suffix ? suffX[index] : Integer.MIN_VALUE));
        finalScores.add(Pair.of(TraceBackState.SUFFIX_Y, noPenaltySeq2Suffix ? suffY[index] : Integer.MIN_VALUE));

        Pair<TraceBackState, Long> bestFinal = Collections.max(finalScores, Comparator.comparingLong(x -> x.getRight()));
        TraceBackState state = bestFinal.getLeft();
        while(r != 0 || c != 0)
        {
            index = (seq2.size() + 1) * r + c;
            if(state == TraceBackState.M)
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
            else if(state == TraceBackState.X)
            {
                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = XTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == TraceBackState.Y)
            {
                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = YTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                c--;
            }
            else if(state == TraceBackState.M2)
            {
                alignment.add(Pair.of(seq1.get(r - 1), seq2.get(c - 1)));
                state = M2Trace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
                c--;
            }
            else if(state == TraceBackState.X2)
            {
                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = X2Trace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == TraceBackState.Y2)
            {
                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = Y2Trace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                c--;
            }
            else if(state == TraceBackState.SUFFIX_X)
            {
                if(!noPenaltySeq1Suffix)
                {
                    throw new RuntimeException("SUFFIX_X is an invalid traceback state");
                }

                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = suffXTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == TraceBackState.SUFFIX_Y)
            {
                if(!noPenaltySeq2Suffix)
                {
                    throw new RuntimeException("SUFFIX_Y is an invalid traceback state");
                }

                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = suffYTrace[index];
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
}
