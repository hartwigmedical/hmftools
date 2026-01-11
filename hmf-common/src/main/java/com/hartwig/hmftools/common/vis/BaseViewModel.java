package com.hartwig.hmftools.common.vis;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.Nullable;

public class BaseViewModel
{
    private static final int MISSING_BASEQ = 0;

    public final boolean IsOverlapped;

    private final Character mCharBase;
    private final SpecialBase mSpecialBase;
    private final int mBaseQ;

    private boolean mIsSoftClip;
    private List<Character> mRightInsertBases;
    private List<Integer> mRightInsertBaseQs;
    private int mRightInsertCount;

    public BaseViewModel(char base)
    {
        mCharBase = base;
        mSpecialBase = null;
        mBaseQ = MISSING_BASEQ;
        mIsSoftClip = false;
        IsOverlapped = false;

        mRightInsertBases = Lists.newArrayList();
        mRightInsertBaseQs = Lists.newArrayList();
        mRightInsertCount = 0;
    }

    public BaseViewModel(char base, int baseQ)
    {
        this(base, baseQ, false, false);
    }

    public BaseViewModel(char base, int baseQ, boolean isSoftClip, boolean isOverlapped)
    {
        mCharBase = base;
        mSpecialBase = null;
        mBaseQ = baseQ;
        mIsSoftClip = isSoftClip;
        IsOverlapped = isOverlapped;

        mRightInsertBases = Lists.newArrayList();
        mRightInsertBaseQs = Lists.newArrayList();
        mRightInsertCount = 0;
    }

    private BaseViewModel(final SpecialBase base, boolean isOverlapped)
    {
        mCharBase = null;
        mSpecialBase = base;
        mBaseQ = MISSING_BASEQ;
        mIsSoftClip = false;
        IsOverlapped = isOverlapped;

        mRightInsertBases = Lists.newArrayList();
        mRightInsertBaseQs = Lists.newArrayList();
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

    @Nullable
    public List<Character> rightInsertBases() { return mRightInsertBases; }
    @Nullable
    public List<Integer> rightInsertBaseQs() { return mRightInsertBaseQs; }

    public void incRightInsertCount(char base, int baseQ)
    {
        if(mRightInsertBases != null)
        {
            mRightInsertBases.add(base);
            mRightInsertBaseQs.add(baseQ);
        }

        ++mRightInsertCount;
    }

    public void incRightInsertCount(int count)
    {
        mRightInsertBases = null;
        mRightInsertBaseQs = null;
        mRightInsertCount += count;
    }

    public void resetRightInsert()
    {
        mRightInsertBases = Lists.newArrayList();
        mRightInsertBaseQs = Lists.newArrayList();
        mRightInsertCount = 0;
    }

    private boolean hasSpecialBase()
    {
        return mSpecialBase != null;
    }

    public boolean isSoftClip() { return mIsSoftClip; }

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

    public void clearSoftClip()
    {
        mIsSoftClip = false;
    }

    private enum SpecialBase
    {
        MISSING,
        DEL
    }
}
