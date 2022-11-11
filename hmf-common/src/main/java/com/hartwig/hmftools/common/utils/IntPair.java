package com.hartwig.hmftools.common.utils;

// Pair of integers. More efficient than using
// org.apache.commons.lang3.tuple.Pair<Integer, Integer> which uses boxing
public class IntPair
{
    public int left;
    public int right;

    public IntPair(int l, int r)
    {
        left = l;
        right = r;
    }

    public int getLeft() { return left; }

    public int getRight() { return right; }

    public int getKey() {
        return this.getLeft();
    }

    public int getValue() {
        return this.getRight();
    }

    @Override
    public boolean equals(final Object o)
    {
        if (this == o)
        {
            return true;
        }
        if (o == null || getClass() != o.getClass())
        {
            return false;
        }
        final IntPair intPair = (IntPair) o;
        return left == intPair.left && right == intPair.right;
    }

    @Override
    public int hashCode()
    {
        return 31 * Integer.hashCode(left) + Integer.hashCode(right);
    }
}
