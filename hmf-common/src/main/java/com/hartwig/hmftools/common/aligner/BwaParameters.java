package com.hartwig.hmftools.common.aligner;

import java.util.List;

import htsjdk.samtools.CigarElement;

public final class BwaParameters
{
    public static final int BWA_MATCH_SCORE = 1;
    public static final int BWA_MISMATCH_PENALTY = 4;
    public static final int BWA_GAP_OPEN_PENALTY = 6;
    public static final int BWA_GAP_EXTEND_PENALTY = 1;
    public static final int BWA_CLIPPING_PENALTY = 5;

    public static int alignmentScore(final List<CigarElement> cigarElements, final int numMutations)
    {
        int mismatchCount = numMutations;
        int alignmentScore = 0;
        for(CigarElement el : cigarElements)
        {
            if(el.getOperator().isClipping())
            {
                alignmentScore -= BWA_CLIPPING_PENALTY;
                continue;
            }

            boolean isRef = el.getOperator().consumesReferenceBases();
            boolean isRead = el.getOperator().consumesReadBases();
            if(isRead ^ isRef)
            {
                alignmentScore -= BWA_GAP_OPEN_PENALTY + el.getLength() * BWA_GAP_EXTEND_PENALTY;
                mismatchCount -= el.getLength();
                continue;
            }

            alignmentScore += el.getLength() * BWA_MATCH_SCORE;
        }

        alignmentScore -= mismatchCount * (BWA_MATCH_SCORE + BWA_MISMATCH_PENALTY);
        return alignmentScore;
    }
}
