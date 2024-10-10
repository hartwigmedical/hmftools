package com.hartwig.hmftools.bamtools.fastqbiomodalcollapse;

import static java.lang.Math.max;

import static com.hartwig.hmftools.bamtools.fastqbiomodalcollapse.NeedlemanWunschAligner.GapTrace.OPEN_GAP;

import java.util.Collections;
import java.util.List;
import java.util.function.ToIntBiFunction;

import org.apache.commons.compress.utils.Lists;
import org.apache.commons.lang3.tuple.Pair;

public class NeedlemanWunschAligner<T>
{
    private final ToIntBiFunction<T, T> mScoringFn;
    private final int mOpenGapPenalty;
    private final int mExtendGapPenalty;

    public NeedlemanWunschAligner(final ToIntBiFunction<T, T> scoringFn, int openGapPenalty, int extendGapPenalty)
    {
        mScoringFn = scoringFn;
        mOpenGapPenalty = openGapPenalty;
        mExtendGapPenalty = extendGapPenalty;
    }

    private enum MatchTrace
    {
        EXTEND_M,
        CLOSE_X_GAP,
        CLOSE_Y_GAP;
    }

    public enum GapTrace
    {
        OPEN_GAP,
        EXTEND_GAP;
    }

    public enum TraceBackState
    {
        M,
        X,
        Y,
        SUFFIX_X,
        SUFFIX_Y;
    }

    public List<Pair<T, T>> align(final List<T> seq1, final List<T> seq2, boolean noPenaltySeq1Suffix, boolean noPenaltySeq2Suffix)
    {
        // initialize scoring matrices
        long[] M = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] X = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] Y = new long[(seq1.size() + 1) * (seq2.size() + 1)];
        long[] suffX = noPenaltySeq1Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        long[] suffY = noPenaltySeq2Suffix ? new long[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        MatchTrace[] MTrace = new MatchTrace[(seq1.size() + 1) * (seq2.size() + 1)];
        GapTrace[] XTrace = new GapTrace[(seq1.size() + 1) * (seq2.size() + 1)];
        GapTrace[] YTrace = new GapTrace[(seq1.size() + 1) * (seq2.size() + 1)];
        GapTrace[] suffXTrace = noPenaltySeq1Suffix ? new GapTrace[(seq1.size() + 1) * (seq2.size() + 1)] : null;
        GapTrace[] suffYTrace = noPenaltySeq2Suffix ? new GapTrace[(seq1.size() + 1) * (seq2.size() + 1)] : null;

        for(int r = 1; r <= seq1.size(); r++)
        {
            int index = (seq2.size() + 1) * r;

            M[index] = Integer.MIN_VALUE;
            X[index] = -mOpenGapPenalty - (r - 1) * mExtendGapPenalty;
            Y[index] = Integer.MIN_VALUE;

            MTrace[index] = null;
            XTrace[index] = (r == 1) ? OPEN_GAP : GapTrace.EXTEND_GAP;
            YTrace[index] = null;

            if(noPenaltySeq1Suffix)
            {
                suffX[index] = 0;
                suffXTrace[index] = (r == 1) ? OPEN_GAP : GapTrace.EXTEND_GAP;
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

            MTrace[c] = null;
            XTrace[c] = null;
            YTrace[c] = (c == 1) ? OPEN_GAP : GapTrace.EXTEND_GAP;

            if(noPenaltySeq1Suffix)
            {
                suffX[c] = Integer.MIN_VALUE;
                suffXTrace[c] = null;
            }

            if(noPenaltySeq2Suffix)
            {
                suffY[c] = 0;
                suffYTrace[c] = (c == 1) ? OPEN_GAP : GapTrace.EXTEND_GAP;
            }
        }

        M[0] = 0;
        X[0] = Integer.MIN_VALUE;
        Y[0] = Integer.MIN_VALUE;

        MTrace[0] = null;
        XTrace[0] = null;
        YTrace[0] = null;

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

                // M
                long extendMatch = M[diagIndex] + matchScore;
                long closeXGap = X[diagIndex] + matchScore;
                long closeYGap = Y[diagIndex] + matchScore;
                long bestMScore = max(extendMatch, max(closeXGap, closeYGap));
                M[index] = bestMScore;

                if(bestMScore == extendMatch)
                {
                    MTrace[index] = MatchTrace.EXTEND_M;
                }
                else if(bestMScore == closeXGap)
                {
                    MTrace[index] = MatchTrace.CLOSE_X_GAP;
                }
                else
                {
                    MTrace[index] = MatchTrace.CLOSE_Y_GAP;
                }

                // X
                long openGap = M[upIndex] - mOpenGapPenalty;
                long extendGap = X[upIndex] - mExtendGapPenalty;
                if(extendGap >= openGap)
                {
                    X[index] = extendGap;
                    XTrace[index] = GapTrace.EXTEND_GAP;
                }
                else
                {
                    X[index] = openGap;
                    XTrace[index] = OPEN_GAP;
                }

                // Y
                openGap = M[leftIndex] - mOpenGapPenalty;
                extendGap = Y[leftIndex] - mExtendGapPenalty;
                if(extendGap >= openGap)
                {
                    Y[index] = extendGap;
                    YTrace[index] = GapTrace.EXTEND_GAP;
                }
                else
                {
                    Y[index] = openGap;
                    YTrace[index] = OPEN_GAP;
                }

                // suffX
                if(noPenaltySeq1Suffix)
                {
                    openGap = M[upIndex];
                    extendGap = suffX[upIndex];
                    if(extendGap >= openGap)
                    {
                        suffX[index] = extendGap;
                        suffXTrace[index] = GapTrace.EXTEND_GAP;
                    }
                    else
                    {
                        suffX[index] = openGap;
                        suffXTrace[index] = OPEN_GAP;
                    }
                }

                // suffY
                if(noPenaltySeq2Suffix)
                {
                    openGap = M[leftIndex];
                    extendGap = suffY[leftIndex];
                    if(extendGap >= openGap)
                    {
                        suffY[index] = extendGap;
                        suffYTrace[index] = GapTrace.EXTEND_GAP;
                    }
                    else
                    {
                        suffY[index] = openGap;
                        suffYTrace[index] = OPEN_GAP;
                    }
                }
            }
        }

        // trace back
        List<Pair<T, T>> alignment = Lists.newArrayList();
        int r = seq1.size();
        int c = seq2.size();
        int index = (seq2.size() + 1) * r + c;

        List<Long> finalScores = Lists.newArrayList();
        finalScores.add(M[index]);
        finalScores.add(X[index]);
        finalScores.add(Y[index]);
        finalScores.add(noPenaltySeq1Suffix ? suffX[index] : Integer.MIN_VALUE);
        finalScores.add(noPenaltySeq2Suffix ? suffY[index] : Integer.MIN_VALUE);

        long bestFinal = Collections.max(finalScores);
        TraceBackState state = null;
        if(bestFinal == M[index])
        {
            state = TraceBackState.M;
        }
        else if(bestFinal == X[index])
        {
            state = TraceBackState.X;
        }
        else if(bestFinal == Y[index])
        {
            state = TraceBackState.Y;
        }
        else if(bestFinal == (noPenaltySeq1Suffix ? suffX[index] : Integer.MIN_VALUE))
        {
            state = TraceBackState.SUFFIX_X;
        }
        else if(bestFinal == (noPenaltySeq2Suffix ? suffY[index] : Integer.MIN_VALUE))
        {
            state = TraceBackState.SUFFIX_Y;
        }
        else
        {
            throw new RuntimeException("Unreachable");
        }

        while(r != 0 || c != 0)
        {
            index = (seq2.size() + 1) * r + c;
            if(state == TraceBackState.M)
            {
                alignment.add(Pair.of(seq1.get(r - 1), seq2.get(c - 1)));
                MatchTrace nextStep = MTrace[index];
                if(nextStep == null)
                {
                    throw new RuntimeException("Invalid next step");
                }

                switch(nextStep)
                {
                    case CLOSE_X_GAP:
                        state = TraceBackState.X;
                        break;
                    case CLOSE_Y_GAP:
                        state = TraceBackState.Y;
                        break;
                }

                r--;
                c--;
            }
            else if(state == TraceBackState.X)
            {
                alignment.add(Pair.of(seq1.get(r - 1), null));
                GapTrace nextStep = XTrace[index];
                if(nextStep == null)
                {
                    throw new RuntimeException("Invalid next step");
                }

                if(nextStep == OPEN_GAP)
                {
                    state = TraceBackState.M;
                }

                r--;
            }
            else if(state == TraceBackState.Y)
            {
                alignment.add(Pair.of(null, seq2.get(c - 1)));
                GapTrace nextStep = YTrace[index];
                if(nextStep == null)
                {
                    throw new RuntimeException("Invalid next step");
                }

                if(nextStep == OPEN_GAP)
                {
                    state = TraceBackState.M;
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
                GapTrace nextStep = suffXTrace[index];
                if(nextStep == null)
                {
                    throw new RuntimeException("Invalid next step");
                }

                if(nextStep == OPEN_GAP)
                {
                    state = TraceBackState.M;
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
                GapTrace nextStep = suffYTrace[index];
                if(nextStep == null)
                {
                    throw new RuntimeException("Invalid next step");
                }

                if(nextStep == OPEN_GAP)
                {
                    state = TraceBackState.M;
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
