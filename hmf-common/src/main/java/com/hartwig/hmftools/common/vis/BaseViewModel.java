package com.hartwig.hmftools.common.vis;

public class BaseViewModel
{
    private static final int MISSING_BASEQ = 0;

    public final boolean IsSoftClip;
    public final boolean IsOverlapped;

    private final Character mCharBase;
    private final SpecialBase mSpecialBase;
    private final int mBaseQ;

    private int mRightInsertCount;

    private BaseViewModel(final Character base, int baseQ, boolean isSoftClip, boolean isOverlapped, int rightInsertCount, final SpecialBase specialBase)
    {
        mCharBase = base;
        mBaseQ = baseQ;
        IsSoftClip = isSoftClip;
        IsOverlapped = isOverlapped;
        mRightInsertCount = rightInsertCount;
        mSpecialBase = specialBase;
    }

    public BaseViewModel(char base)
    {
        mCharBase = base;
        mSpecialBase = null;
        mBaseQ = MISSING_BASEQ;
        IsSoftClip = false;
        IsOverlapped = false;

        mRightInsertCount = 0;
    }

    public BaseViewModel(char base, int baseQ, boolean isSoftClip, boolean isOverlapped)
    {
        mCharBase = base;
        mSpecialBase = null;
        mBaseQ = baseQ;
        IsSoftClip = isSoftClip;
        IsOverlapped = isOverlapped;

        mRightInsertCount = 0;
    }

    private BaseViewModel(final SpecialBase base, boolean isOverlapped)
    {
        mCharBase = null;
        mSpecialBase = base;
        mBaseQ = MISSING_BASEQ;
        IsSoftClip = false;
        IsOverlapped = isOverlapped;

        mRightInsertCount = 0;
    }

    public static BaseViewModel createMissingBase()
    {
        return new BaseViewModel(SpecialBase.MISSING, false);
    }

    public static BaseViewModel createDelBase()
    {
        return new BaseViewModel(SpecialBase.DEL, false);
    }

    public static BaseViewModel createDelBase(boolean isOverlapped)
    {
        return new BaseViewModel(SpecialBase.DEL, isOverlapped);
    }

    public BaseViewModel clearSoftClip()
    {
        return new BaseViewModel(mCharBase, mBaseQ, false, IsOverlapped, mRightInsertCount, mSpecialBase);
    }

    public char charBase()
    {
        if(!hasCharBase())
        {
            throw new IllegalStateException("BaseViewModel does not have character base.");
        }

        return mCharBase;
    }

    private SpecialBase specialBase()
    {
        if(!hasSpecialBase())
        {
            throw new IllegalStateException("BaseViewModel does not have special base.");
        }

        return mSpecialBase;
    }

    public int baseQ()
    {
        return mBaseQ;
    }

    public int rightInsertCount()
    {
        return mRightInsertCount;
    }

    public void incRightInsertCount()
    {
        ++mRightInsertCount;
    }

    public void incRightInsertCount(int count)
    {
        mRightInsertCount += count;
    }

    private boolean hasSpecialBase()
    {
        return mSpecialBase != null;
    }

    public boolean hasCharBase()
    {
        return mCharBase != null;
    }

    public boolean isMissing()
    {
        return hasSpecialBase() && specialBase() == SpecialBase.MISSING;
    }

    public boolean isDel()
    {
        return hasSpecialBase() && specialBase() == SpecialBase.DEL;
    }

    private enum SpecialBase
    {
        MISSING,
        DEL
    }
}
