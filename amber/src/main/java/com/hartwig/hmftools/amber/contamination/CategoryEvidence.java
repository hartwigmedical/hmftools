package com.hartwig.hmftools.amber.contamination;

import java.util.Objects;

import com.google.common.base.Preconditions;

import org.jetbrains.annotations.NotNull;

public class CategoryEvidence<T extends Comparable<T>> implements Comparable<CategoryEvidence<T>>
{
    private final T Category;
    private int TotalPoints = 0;
    private int EvidencePoints = 0;

    public CategoryEvidence(final T Category)
    {
        Preconditions.checkNotNull(Category);
        this.Category = Category;
    }

    public T category()
    {
        return Category;
    }

    public void register(boolean isEvidence)
    {
        TotalPoints++;
        if(isEvidence)
        {
            EvidencePoints++;
        }
    }

    public void set(int totalPoints, int evidencePoints)
    {
        TotalPoints = totalPoints;
        EvidencePoints = evidencePoints;
    }

    public int totalPoints()
    {
        return TotalPoints;
    }

    public int evidencePoints()
    {
        return EvidencePoints;
    }

    @Override
    public int compareTo(@NotNull final CategoryEvidence<T> o)
    {
        int result = Double.compare(ratio(), o.ratio());
        if(result == 0)
        {
            result = Category.compareTo(o.Category);
        }
        return result;
    }

    public double ratio()
    {
        if(TotalPoints == 0)
        {
            return Double.MAX_VALUE;
        }
        return (double) EvidencePoints / (double) TotalPoints;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final CategoryEvidence<?> that = (CategoryEvidence<?>) o;
        return Objects.equals(Category, that.Category);
    }

    @Override
    public int hashCode()
    {
        return Objects.hashCode(Category);
    }

    @Override
    public String toString()
    {
        return String.format("%s: %d/%d", Category, EvidencePoints, TotalPoints);
    }
}
