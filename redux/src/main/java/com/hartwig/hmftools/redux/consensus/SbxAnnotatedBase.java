package com.hartwig.hmftools.redux.consensus;

import static java.lang.String.format;

import static htsjdk.samtools.CigarOperator.S;

import htsjdk.samtools.CigarOperator;

public class SbxAnnotatedBase
{
    public final int ReadIndex;
    public final int RefPos;
    public final CigarOperator Op;
    public final byte ReadBase;
    public final boolean IsDuplexIndel;

    private byte mQual;
    private boolean mDeleted;
    private boolean mInSoftClip; // true if derived from a soft-clip but has been replaced by the aligned supplementary data

    protected final static byte INVALID_BASE = -1;

    public SbxAnnotatedBase(int readIndex, int refPos, final CigarOperator op, byte readBase, byte qual, boolean isDuplexIndel)
    {
        ReadIndex = readIndex;
        RefPos = refPos;
        Op = op;
        ReadBase = readBase;
        IsDuplexIndel = isDuplexIndel;

        mQual = qual;
        mDeleted = false;
        mInSoftClip = false;
    }

    public boolean isReadBase() { return Op.consumesReadBases(); }
    public boolean isRefBase()
    {
        return Op.consumesReferenceBases();
    }
    public void deleteBase()
    {
        mDeleted = true;
    }
    public boolean deleted()
    {
        return mDeleted;
    }
    public byte qual()
    {
        return mQual;
    }

    public void setInSoftClip() { mInSoftClip = true; }

    public boolean inSoftClip() { return Op == S || mInSoftClip; }
    public CigarOperator originalOperator() { return mInSoftClip ? S : Op; }

    public boolean setQual(byte qual)
    {
        if(mQual == qual)
            return false;

        mQual = qual;
        return true;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;

        if(!(o instanceof SbxAnnotatedBase))
            return false;

        final SbxAnnotatedBase that = (SbxAnnotatedBase) o;
        return ReadIndex == that.ReadIndex && RefPos == that.RefPos && ReadBase == that.ReadBase && IsDuplexIndel == that.IsDuplexIndel
                && mQual == that.mQual && mDeleted == that.mDeleted && Op == that.Op;
    }

    @Override
    public int hashCode()
    {
        int hash = ReadIndex;
        hash = 31 * hash + RefPos;
        hash = 31 * hash + Op.hashCode();
        hash = 31 * hash + (int) ReadBase;
        hash = 31 * hash + (IsDuplexIndel ? 1 : 0);
        hash = 31 * hash + (int) mQual;
        hash = 31 * hash + (mDeleted ? 1 : 0);

        return hash;
    }

    @Override
    public String toString()
    {
        return format("%d:%d %c@%d cigar(%s) duplexIndel(%s) deleted(%s) softClip(%s)",
                ReadIndex, RefPos, (char)ReadBase, mQual, Op, IsDuplexIndel, mDeleted, mInSoftClip);
    }
}
