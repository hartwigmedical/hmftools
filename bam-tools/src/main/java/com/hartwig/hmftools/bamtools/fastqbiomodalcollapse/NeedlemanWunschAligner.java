package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.ToIntBiFunction;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.Nullable;

public class NeedlemanWunschAligner<T>
{

    public NeedlemanWunschAligner()
    {
    }

    public enum TraceBackState
    {
        M,
        X,
        Y,
        PREFIX_X,
        PREFIX_Y,
        SUFFIX_X,
        SUFFIX_Y;
    }

    public List<Pair<T, T>> align(final List<T> seq1, final List<T> seq2, final ToIntBiFunction<T, T> scoringFn, int openGapPenalty,
            int extendGapPenalty, boolean noPenaltySeq1Prefix, boolean noPenaltySeq2Prefix, boolean noPenaltySeq1Suffix,
            boolean noPenaltySeq2Suffix)
    {
        // initialize scoring matrices
        long[] M = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] X = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] Y = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] preX = noPenaltySeq1Prefix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        long[] preY = noPenaltySeq2Prefix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        long[] suffX = noPenaltySeq1Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        long[] suffY = noPenaltySeq2Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        TraceBackState[] MTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] XTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] YTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
        TraceBackState[] preXTrace = noPenaltySeq1Prefix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        TraceBackState[] preYTrace = noPenaltySeq2Prefix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        TraceBackState[] suffXTrace = noPenaltySeq1Suffix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        TraceBackState[] suffYTrace = noPenaltySeq2Suffix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        for(int r = 1; r <= seq1.size(); r++)
        {
            int index = (seq2.size() + 1) * r;

            M[index] = Integer.MIN_VALUE;
            X[index] = -openGapPenalty - (r - 1) * extendGapPenalty;
            Y[index] = Integer.MIN_VALUE;

            MTrace[index] = null;
            XTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.X;
            YTrace[index] = null;

            if(noPenaltySeq1Prefix)
            {
                preX[index] = 0;
                preXTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.PREFIX_X;
            }

            if(noPenaltySeq2Prefix)
            {
                preY[index] = Integer.MIN_VALUE;
                preYTrace[index] = null;
            }

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
            Y[c] = -openGapPenalty - (c - 1) * extendGapPenalty;

            MTrace[c] = null;
            XTrace[c] = null;
            YTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.Y;

            if(noPenaltySeq1Prefix)
            {
                preX[c] = Integer.MIN_VALUE;
                preXTrace[c] = null;
            }

            if(noPenaltySeq2Prefix)
            {
                preY[c] = 0;
                preYTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.PREFIX_Y;
            }

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

        MTrace[0] = null;
        XTrace[0] = null;
        YTrace[0] = null;

        if(noPenaltySeq1Prefix)
        {
            preX[0] = Integer.MIN_VALUE;
            preXTrace[0] = null;
        }

        if(noPenaltySeq2Prefix)
        {
            preY[0] = Integer.MIN_VALUE;
            preYTrace[0] = null;
        }

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
        List<Pair<TraceBackState, Long>> edgeScores = Lists.newArrayList();
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

                // M
                int matchScore = scoringFn.applyAsInt(base1, base2);
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[diagIndex] + matchScore));
                edgeScores.add(Pair.of(TraceBackState.X, X[diagIndex] + matchScore));
                edgeScores.add(Pair.of(TraceBackState.Y, Y[diagIndex] + matchScore));
                edgeScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltySeq1Prefix ? preX[diagIndex] + matchScore : Integer.MIN_VALUE));
                edgeScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltySeq2Prefix ? preY[diagIndex] + matchScore : Integer.MIN_VALUE));
                Pair<TraceBackState, Long> bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                M[index] = bestTransition.getRight();
                MTrace[index] = bestTransition.getLeft();

                // X
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[upIndex] - openGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.X, X[upIndex] - extendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                X[index] = bestTransition.getRight();
                XTrace[index] = bestTransition.getLeft();

                // Y
                edgeScores.clear();
                edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex] - openGapPenalty));
                edgeScores.add(Pair.of(TraceBackState.Y, Y[leftIndex] - extendGapPenalty));
                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                Y[index] = bestTransition.getRight();
                YTrace[index] = bestTransition.getLeft();

                if(noPenaltySeq1Prefix)
                {
                    // preX
                    preX[index] = preX[upIndex];
                    preXTrace[index] = TraceBackState.PREFIX_X;
                }

                if(noPenaltySeq2Prefix)
                {
                    // preX
                    preY[index] = preY[leftIndex];
                    preYTrace[index] = TraceBackState.PREFIX_Y;
                }

                if(noPenaltySeq1Suffix)
                {
                    // suffX
                    edgeScores.clear();
                    edgeScores.add(Pair.of(TraceBackState.M, M[upIndex]));
                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_X, suffX[upIndex]));
                    edgeScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltySeq2Prefix ? preY[upIndex] : Integer.MIN_VALUE));
                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
                    suffX[index] = bestTransition.getRight();
                    suffXTrace[index] = bestTransition.getLeft();
                }

                if(noPenaltySeq2Suffix)
                {
                    // suffY
                    edgeScores.clear();
                    edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex]));
                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_Y, suffY[leftIndex]));
                    edgeScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltySeq1Prefix ? preX[leftIndex] : Integer.MIN_VALUE));
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
        finalScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltySeq1Prefix ? preX[index] : Integer.MIN_VALUE));
        finalScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltySeq2Prefix ? preY[index] : Integer.MIN_VALUE));
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
            else if(state == TraceBackState.PREFIX_X)
            {
                if(!noPenaltySeq1Prefix)
                {
                    throw new RuntimeException("PREFIX_X is an invalid traceback state");
                }

                alignment.add(Pair.of(seq1.get(r - 1), null));
                state = preXTrace[index];
                if(state == null)
                {
                    throw new RuntimeException("Invalid traceback transition");
                }

                r--;
            }
            else if(state == TraceBackState.PREFIX_Y)
            {
                if(!noPenaltySeq2Prefix)
                {
                    throw new RuntimeException("PREFIX_Y is an invalid traceback state");
                }

                alignment.add(Pair.of(null, seq2.get(c - 1)));
                state = preYTrace[index];
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

//    public List<Pair<T, T>> align(final List<T> seq1, final List<T> seq2, final ToIntBiFunction<T, T> scoringFn, int openGapPenalty,
//            int extendGapPenalty, boolean noPenaltyPrefix, final SuffixPenaltyOption suffixPenaltyOption)
//    {
//        return align2(seq1, seq2, scoringFn, null, openGapPenalty, extendGapPenalty, noPenaltyPrefix, suffixPenaltyOption);
//    }
//
//    public List<Pair<T, T>> align2(final List<T> seq1, final List<T> seq2, final ToIntBiFunction<T, T> scoringFn,
//            @Nullable final ToIntBiFunction<T, T> scoringFn2, int openGapPenalty, int extendGapPenalty, boolean noPenaltyPrefix,
//            final SuffixPenaltyOption suffixPenaltyOption)
//    {
//        // initialize scoring matrices
//        long[] M = new long[(seq1.size() + 1) * (seq2.size() + 1)];
//        long[] X = new long[(seq1.size() + 1) * (seq2.size() + 1)];
//        long[] Y = new long[(seq1.size() + 1) * (seq2.size() + 1)];
//        long[] M2 = scoringFn2 != null ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] X2 = scoringFn2 != null ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] Y2 = scoringFn2 != null ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] preX = noPenaltyPrefix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] preY = noPenaltyPrefix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] suffX = suffixPenaltyOption != PENALIZE ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        long[] suffY = suffixPenaltyOption != PENALIZE ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//
//        TraceBackState[] MTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
//        TraceBackState[] XTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
//        TraceBackState[] YTrace = new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)];
//        TraceBackState[] M2Trace = scoringFn2 != null ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] X2Trace = scoringFn2 != null ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] Y2Trace = scoringFn2 != null ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] preXTrace = noPenaltyPrefix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] preYTrace = noPenaltyPrefix ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] suffXTrace = suffixPenaltyOption != PENALIZE ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//        TraceBackState[] suffYTrace = suffixPenaltyOption != PENALIZE ? new TraceBackState[(seq1.size() + 1) * (seq2.size() + 1)] : null;
//
//        for(int r = 1; r <= seq1.size(); r++)
//        {
//            int index = (seq2.size() + 1) * r;
//
//            M[index] = Integer.MIN_VALUE;
//            X[index] = -openGapPenalty - (r - 1) * extendGapPenalty;
//            Y[index] = Integer.MIN_VALUE;
//
//            MTrace[index] = null;
//            XTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.X;
//            YTrace[index] = null;
//
//            if(scoringFn2 != null)
//            {
//                M2[index] = Integer.MIN_VALUE;
//                X2[index] = Integer.MIN_VALUE;
//                Y2[index] = Integer.MIN_VALUE;
//
//                M2Trace[index] = null;
//                X2Trace[index] = null;
//                Y2Trace[index] = null;
//            }
//
//            if(noPenaltyPrefix)
//            {
//                preX[index] = 0;
//                preXTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.PREFIX_X;
//
//                preY[index] = Integer.MIN_VALUE;
//                preYTrace[index] = null;
//            }
//
//            if(suffixPenaltyOption != PENALIZE)
//            {
//                suffX[index] = 0;
//                suffXTrace[index] = (r == 1) ? TraceBackState.M : TraceBackState.SUFFIX_X;
//
//                suffY[index] = Integer.MIN_VALUE;
//                suffYTrace[index] = null;
//            }
//        }
//
//        for(int c = 1; c <= seq2.size(); c++)
//        {
//            M[c] = Integer.MIN_VALUE;
//            X[c] = Integer.MIN_VALUE;
//            Y[c] = -openGapPenalty - (c - 1) * extendGapPenalty;
//
//            MTrace[c] = null;
//            XTrace[c] = null;
//            YTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.Y;
//
//            if(scoringFn2 != null)
//            {
//                M2[c] = Integer.MIN_VALUE;
//                X2[c] = Integer.MIN_VALUE;
//                Y2[c] = Integer.MIN_VALUE;
//
//                M2Trace[c] = null;
//                X2Trace[c] = null;
//                Y2Trace[c] = null;
//            }
//
//            if(noPenaltyPrefix)
//            {
//                preX[c] = Integer.MIN_VALUE;
//                preXTrace[c] = null;
//
//                preY[c] = 0;
//                preYTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.PREFIX_Y;
//            }
//
//            if(suffixPenaltyOption != PENALIZE)
//            {
//                suffX[c] = Integer.MIN_VALUE;
//                suffXTrace[c] = null;
//
//                suffY[c] = 0;
//                suffYTrace[c] = (c == 1) ? TraceBackState.M : TraceBackState.SUFFIX_Y;
//            }
//        }
//
//        M[0] = 0;
//        X[0] = Integer.MIN_VALUE;
//        Y[0] = Integer.MIN_VALUE;
//
//        MTrace[0] = null;
//        XTrace[0] = null;
//        YTrace[0] = null;
//
//        if(scoringFn2 != null)
//        {
//            M2[0] = Integer.MIN_VALUE;
//            X2[0] = Integer.MIN_VALUE;
//            Y2[0] = Integer.MIN_VALUE;
//
//            M2Trace[0] = null;
//            X2Trace[0] = null;
//            Y2Trace[0] = null;
//        }
//
//        if(noPenaltyPrefix)
//        {
//            preX[0] = Integer.MIN_VALUE;
//            preXTrace[0] = null;
//
//            preY[0] = Integer.MIN_VALUE;
//            preYTrace[0] = null;
//        }
//
//        if(suffixPenaltyOption != PENALIZE)
//        {
//            suffX[0] = Integer.MIN_VALUE;
//            suffXTrace[0] = null;
//
//            suffY[0] = Integer.MIN_VALUE;
//            suffYTrace[0] = null;
//        }
//
//        // fill out
//        List<Pair<TraceBackState, Long>> edgeScores = Lists.newArrayList();
//        for(int r = 1; r <= seq1.size(); r++)
//        {
//            T base1 = seq1.get(r - 1);
//            for(int c = 1; c <= seq2.size(); c++)
//            {
//                T base2 = seq2.get(c - 1);
//                int index = (seq2.size() + 1) * r + c;
//                int upIndex = (seq2.size() + 1) * (r - 1) + c;
//                int leftIndex = (seq2.size() + 1) * r + c - 1;
//                int diagIndex = (seq2.size() + 1) * (r - 1) + c - 1;
//
//                // M
//                int matchScore = scoringFn.applyAsInt(base1, base2);
//                edgeScores.clear();
//                edgeScores.add(Pair.of(TraceBackState.M, M[diagIndex] + matchScore));
//                edgeScores.add(Pair.of(TraceBackState.X, X[diagIndex] + matchScore));
//                edgeScores.add(Pair.of(TraceBackState.Y, Y[diagIndex] + matchScore));
//                edgeScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltyPrefix ? preX[diagIndex] + matchScore : Integer.MIN_VALUE));
//                edgeScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltyPrefix ? preY[diagIndex] + matchScore : Integer.MIN_VALUE));
//                Pair<TraceBackState, Long> bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                M[index] = bestTransition.getRight();
//                MTrace[index] = bestTransition.getLeft();
//
//                // X
//                edgeScores.clear();
//                edgeScores.add(Pair.of(TraceBackState.M, M[upIndex] - openGapPenalty));
//                edgeScores.add(Pair.of(TraceBackState.X, X[upIndex] - extendGapPenalty));
//                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                X[index] = bestTransition.getRight();
//                XTrace[index] = bestTransition.getLeft();
//
//                // Y
//                edgeScores.clear();
//                edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex] - openGapPenalty));
//                edgeScores.add(Pair.of(TraceBackState.Y, Y[leftIndex] - extendGapPenalty));
//                bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                Y[index] = bestTransition.getRight();
//                YTrace[index] = bestTransition.getLeft();
//
//                if(scoringFn2 != null)
//                {
//                    // M2
//                    int matchScore2 = scoringFn2.applyAsInt(base1, base2);
//                    edgeScores.clear();
//                    edgeScores.add(Pair.of(TraceBackState.M, M[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.X, X[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.Y, Y[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.M2, M2[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.X2, X2[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.Y2, Y2[diagIndex] + matchScore2));
//                    edgeScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltyPrefix ? preX[diagIndex] + matchScore2 : Integer.MIN_VALUE));
//                    edgeScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltyPrefix ? preY[diagIndex] + matchScore2 : Integer.MIN_VALUE));
//                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                    M2[index] = bestTransition.getRight();
//                    M2Trace[index] = bestTransition.getLeft();
//
//                    // X2
//                    edgeScores.clear();
//                    edgeScores.add(Pair.of(TraceBackState.M2, M2[upIndex] - openGapPenalty));
//                    edgeScores.add(Pair.of(TraceBackState.X2, X2[upIndex] - extendGapPenalty));
//                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                    X2[index] = bestTransition.getRight();
//                    X2Trace[index] = bestTransition.getLeft();
//
//                    // Y2
//                    edgeScores.clear();
//                    edgeScores.add(Pair.of(TraceBackState.M2, M2[leftIndex] - openGapPenalty));
//                    edgeScores.add(Pair.of(TraceBackState.Y2, Y2[leftIndex] - extendGapPenalty));
//                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                    Y2[index] = bestTransition.getRight();
//                    Y2Trace[index] = bestTransition.getLeft();
//                }
//
//                if(noPenaltyPrefix)
//                {
//                    // preX
//                    preX[index] = preX[upIndex];
//                    preXTrace[index] = TraceBackState.PREFIX_X;
//
//                    // preX
//                    preY[index] = preY[leftIndex];
//                    preYTrace[index] = TraceBackState.PREFIX_Y;
//                }
//
//                if(suffixPenaltyOption != PENALIZE)
//                {
//                    // suffX
//                    edgeScores.clear();
//                    edgeScores.add(Pair.of(TraceBackState.M, M[upIndex]));
//                    edgeScores.add(Pair.of(TraceBackState.M2, scoringFn2 != null ? M2[upIndex] : Integer.MIN_VALUE));
//                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_X, suffX[upIndex]));
//                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_Y,
//                            suffixPenaltyOption == NO_PENALTY_DUAL ? suffY[upIndex] : Integer.MIN_VALUE));
//                    edgeScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltyPrefix ? preY[upIndex] : Integer.MIN_VALUE));
//                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                    suffX[index] = bestTransition.getRight();
//                    suffXTrace[index] = bestTransition.getLeft();
//
//                    // suffY
//                    edgeScores.clear();
//                    edgeScores.add(Pair.of(TraceBackState.M, M[leftIndex]));
//                    edgeScores.add(Pair.of(TraceBackState.M2, scoringFn2 != null ? M2[leftIndex] : Integer.MIN_VALUE));
//                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_Y, suffY[leftIndex]));
//                    edgeScores.add(Pair.of(TraceBackState.SUFFIX_X,
//                            suffixPenaltyOption == NO_PENALTY_DUAL ? suffX[leftIndex] : Integer.MIN_VALUE));
//                    edgeScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltyPrefix ? preX[leftIndex] : Integer.MIN_VALUE));
//                    bestTransition = Collections.max(edgeScores, Comparator.comparingLong(x -> x.getRight()));
//                    suffY[index] = bestTransition.getRight();
//                    suffYTrace[index] = bestTransition.getLeft();
//                }
//            }
//        }
//
//        // trace back
//        List<Pair<T, T>> alignment = Lists.newArrayList();
//        int r = seq1.size();
//        int c = seq2.size();
//        int index = (seq2.size() + 1) * r + c;
//
//        List<Pair<TraceBackState, Long>> finalScores = Lists.newArrayList();
//        finalScores.add(Pair.of(TraceBackState.M, M[index]));
//        finalScores.add(Pair.of(TraceBackState.X, X[index]));
//        finalScores.add(Pair.of(TraceBackState.Y, Y[index]));
//        finalScores.add(Pair.of(TraceBackState.M2, scoringFn2 != null ? M2[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.X2, scoringFn2 != null ? X2[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.Y2, scoringFn2 != null ? Y2[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.PREFIX_X, noPenaltyPrefix ? preX[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.PREFIX_Y, noPenaltyPrefix ? preY[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.SUFFIX_X, suffixPenaltyOption != PENALIZE ? suffX[index] : Integer.MIN_VALUE));
//        finalScores.add(Pair.of(TraceBackState.SUFFIX_Y, suffixPenaltyOption != PENALIZE ? suffY[index] : Integer.MIN_VALUE));
//
//        Pair<TraceBackState, Long> bestFinal = Collections.max(finalScores, Comparator.comparingLong(x -> x.getRight()));
//        TraceBackState state = bestFinal.getLeft();
//        while(r != 0 || c != 0)
//        {
//            index = (seq2.size() + 1) * r + c;
//            if(state == TraceBackState.M)
//            {
//                alignment.add(Pair.of(seq1.get(r - 1), seq2.get(c - 1)));
//                state = MTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//                c--;
//            }
//            else if(state == TraceBackState.X)
//            {
//                alignment.add(Pair.of(seq1.get(r - 1), null));
//                state = XTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//            }
//            else if(state == TraceBackState.Y)
//            {
//                alignment.add(Pair.of(null, seq2.get(c - 1)));
//                state = YTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                c--;
//            }
//            else if(state == TraceBackState.M2)
//            {
//                if(scoringFn2 == null)
//                {
//                    throw new RuntimeException("M2 is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(seq1.get(r - 1), seq2.get(c - 1)));
//                state = M2Trace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//                c--;
//            }
//            else if(state == TraceBackState.X2)
//            {
//                if(scoringFn2 == null)
//                {
//                    throw new RuntimeException("X2 is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(seq1.get(r - 1), null));
//                state = X2Trace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//            }
//            else if(state == TraceBackState.Y2)
//            {
//                if(scoringFn2 == null)
//                {
//                    throw new RuntimeException("Y2 is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(null, seq2.get(c - 1)));
//                state = Y2Trace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                c--;
//            }
//            else if(state == TraceBackState.PREFIX_X)
//            {
//                if(!noPenaltyPrefix)
//                {
//                    throw new RuntimeException("PREFIX_X is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(seq1.get(r - 1), null));
//                state = preXTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//            }
//            else if(state == TraceBackState.PREFIX_Y)
//            {
//                if(!noPenaltyPrefix)
//                {
//                    throw new RuntimeException("PREFIX_Y is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(null, seq2.get(c - 1)));
//                state = preYTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                c--;
//            }
//            else if(state == TraceBackState.SUFFIX_X)
//            {
//                if(suffixPenaltyOption == PENALIZE)
//                {
//                    throw new RuntimeException("SUFFIX_X is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(seq1.get(r - 1), null));
//                state = suffXTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                r--;
//            }
//            else if(state == TraceBackState.SUFFIX_Y)
//            {
//                if(suffixPenaltyOption == PENALIZE)
//                {
//                    throw new RuntimeException("SUFFIX_Y is an invalid traceback state");
//                }
//
//                alignment.add(Pair.of(null, seq2.get(c - 1)));
//                state = suffYTrace[index];
//                if(state == null)
//                {
//                    throw new RuntimeException("Invalid traceback transition");
//                }
//
//                c--;
//            }
//            else
//            {
//                throw new RuntimeException("Unreachable");
//            }
//        }
//
//        Collections.reverse(alignment);
//        return alignment;
//    }
}
