package com.hartwig.hmftools.svassembly.mermaid;

import org.jetbrains.annotations.Nullable;

class Link<T>
{
    public T Left, Right;
    @Nullable
    public String Text;

    public Link(final T left, final T right, @Nullable final String text)
    {
        Left = left;
        Right = right;
        Text = text;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(this == o)
            return true;
        if(o == null || getClass() != o.getClass())
            return false;

        final Link<?> link = (Link<?>) o;
        return Left == link.Left && Right == link.Right;
    }

    @Override
    public int hashCode()
    {
        return System.identityHashCode(Left) ^ System.identityHashCode(Right);
    }
}
