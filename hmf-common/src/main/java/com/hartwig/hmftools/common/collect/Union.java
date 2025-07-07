package com.hartwig.hmftools.common.collect;

import org.jetbrains.annotations.Nullable;

public class Union<U, V>
{
    private final U mLeft;
    private final V mRight;
    private boolean mIsLeft;

    private Union(final U left, final V right, boolean isLeft)
    {
        mLeft = left;
        mRight = right;
        mIsLeft = isLeft;
    }

    public static <U, V> Union<U, V> createLeft(final U value)
    {
        return new Union<>(value, null, true);
    }

    public static <U, V> Union<U, V> createRight(final V value)
    {
        return new Union<>(null, value, false);
    }

    public boolean hasLeft() { return mIsLeft; }
    public boolean hasRight() { return !mIsLeft; }
    @Nullable public U left() { return mLeft; }
    @Nullable public V right() { return mRight; }
}
