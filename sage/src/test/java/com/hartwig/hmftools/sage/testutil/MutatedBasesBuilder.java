package com.hartwig.hmftools.sage.testutil;

import static htsjdk.samtools.CigarOperator.EQ;
import static htsjdk.samtools.CigarOperator.I;
import static htsjdk.samtools.CigarOperator.X;

import java.util.List;

import com.google.common.collect.Lists;

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
        {
            mMutatedBases[refPos - 1 + i] = null;
        }
    }

    public void insertBases(int refPos, final String bases)
    {
        mMutatedBases[refPos - 1].Insert.append(bases);
    }

    public MutatedBases build()
    {
        List<MutatedBases.MutatedBase> mutatedBases = Lists.newArrayList();
        for(int i = 0; i < mMutatedBases.length; i++)
        {
            if(mMutatedBases[i] == null)
            {
                continue;
            }

            int refPos = i + 1;
            char base = mMutatedBases[i].Base;
            CigarOperator cigarOp = mRefBases.charAt(refPos - 1) == base ? EQ : X;
            mutatedBases.add(new MutatedBases.MutatedBase(refPos, base, cigarOp));

            String insert = mMutatedBases[i].Insert.toString();
            for(int j = 0; j < insert.length(); ++j)
            {
                mutatedBases.add(new MutatedBases.MutatedBase(refPos, insert.charAt(j), I));
            }
        }

        return new MutatedBases(mRefBases, mutatedBases);
    }
}
