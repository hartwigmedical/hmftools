package com.hartwig.hmftools.amber.purity;

import java.util.Objects;

import com.google.common.base.Preconditions;

public class CategoryEvidence<T extends Comparable<T>> implements Comparable<CategoryEvidence<T>>
{
    private final T mCategory;
    private int mTotalPoints = 0;
    private int mEvidencePoints = 0;

    public CategoryEvidence(final T Category)
    {
        Preconditions.checkNotNull(Category);
        this.mCategory = Category;
    }

    public void register(boolean isEvidence)
    {
        mTotalPoints++;
        if(isEvidence)
        {
            mEvidencePoints++;
        }
    }

    public void set(int totalPoints, int evidencePoints)
    {
        mTotalPoints = totalPoints;
        mEvidencePoints = evidencePoints;
    }

    public int totalPoints()
    {
        return mTotalPoints;
    }

    public int evidencePoints()
    {
        return mEvidencePoints;
    }

    @Override
    public int compareTo(final CategoryEvidence<T> o)
    {
        int result = Double.compare(ratio(), o.ratio());
        if(result == 0)
        {
            result = mCategory.compareTo(o.mCategory);
        }
        return result;
    }

    public double ratio()
    {
        if(mTotalPoints == 0)
        {
            return Double.MAX_VALUE;
        }
        return (double) mEvidencePoints / (double) mTotalPoints;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CategoryEvidence<?> that = (CategoryEvidence<?>) o;
        return Objects.equals(mCategory, that.mCategory);
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(mCategory);
    }

    @Override
    public String toString()
    {
        return String.format("%s: %d/%d", mCategory, mEvidencePoints, mTotalPoints);
    }
}
