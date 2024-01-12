package com.hartwig.hmftools.sage.testutil;

import static java.lang.Math.max;
import static java.lang.Math.min;

import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.BaseRegion;

import htsjdk.samtools.CigarOperator;

public class MutatedBasesBuilder
{
    private static class BaseWithInsert
    {
        public char Base;
        public StringBuilder Insert;

        public BaseWithInsert(char base)
        {
            Base = base;
            Insert = new StringBuilder();
        }
    }

    private final String mRefBases;
    private final BaseWithInsert[] mMutatedBases;

    public MutatedBasesBuilder(final String refBases)
    {
        mRefBases = refBases;

        mMutatedBases = new BaseWithInsert[refBases.length()];
        for(int i = 0; i < refBases.length(); ++i)
        {
            char refBase = refBases.charAt(i);
            mMutatedBases[i] = new BaseWithInsert(refBase);
        }
    }

    public void mutateBase(int refPos, char base)
    {
        mMutatedBases[refPos - 1].Base = base;
    }

    public void delBases(int refPos, int length)
    {
        for(int i = 0; i < length; i++)
            mMutatedBases[refPos - 1 + i] = null;
    }

    public void insertBases(int refPos, final String bases)
    {
        mMutatedBases[refPos - 1].Insert.append(bases);
    }

    public MutatedBases build()
    {
        List<MutatedBases.MutatedBase> mutatedBases = Lists.newArrayList();
        int minMutationBound = Integer.MAX_VALUE;
        int maxMutationBound = Integer.MIN_VALUE;
        for(int i = 0; i < mMutatedBases.length; i++)
        {
            int refPos = i + 1;
            if(mMutatedBases[i] == null)
            {
                minMutationBound = min(minMutationBound, refPos);
                maxMutationBound = max(maxMutationBound, refPos);
                continue;
            }

            char base = mMutatedBases[i].Base;
            CigarOperator cigarOp = mRefBases.charAt(refPos - 1) == base ? EQ : X;
            if(mRefBases.charAt(refPos - 1) != base)
            {
                minMutationBound = min(minMutationBound, refPos);
                maxMutationBound = max(maxMutationBound, refPos);
            }

            mutatedBases.add(new MutatedBases.MutatedBase(refPos, base, cigarOp));

            String insert = mMutatedBases[i].Insert.toString();
            for(int j = 0; j < insert.length(); ++j)
                mutatedBases.add(new MutatedBases.MutatedBase(refPos, insert.charAt(j), I));

            if(!insert.isEmpty())
            {
                minMutationBound = min(minMutationBound, refPos);
                maxMutationBound = max(maxMutationBound, refPos + 1);
            }
        }

        BaseRegion mutationBounds = minMutationBound <= maxMutationBound ? new BaseRegion(minMutationBound, maxMutationBound) : null;
        return new MutatedBases(mRefBases, mutatedBases, mutationBounds);
    }
}
