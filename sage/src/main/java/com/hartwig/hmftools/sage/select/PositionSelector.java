package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.region.BasePosition;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class PositionSelector<P extends BasePosition>
{
    @NotNull
    private final List<P> mPositions;
    private int mIndex = 0;

    public PositionSelector(@NotNull final List<P> positions)
    {
        this.mPositions = positions;
    }

    @NotNull
    public Optional<P> select(long position)
    {
        if(mPositions.isEmpty())
        {
            return Optional.empty();
        }

        int currentCompare = Long.compare(current().Position, position);
        while(currentCompare >= 0 && mIndex > 0)
        {
            mIndex--;
            currentCompare = Long.compare(current().Position, position);
        }

        while(currentCompare < 0 && mIndex < mPositions.size() - 1)
        {
            mIndex++;
            currentCompare = Long.compare(current().Position, position);
        }

        return currentCompare == 0 ? Optional.of(current()) : Optional.empty();
    }

    public void select(final int start, final int end, final Consumer<P> handler)
    {
        if(mPositions.isEmpty())
        {
            return;
        }

        // Line up index to start
        select(start);

        while(mIndex < mPositions.size() && inRegion(current(), start, end))
        {
            handler.accept(current());
            mIndex++;
        }

        mIndex = Math.min(mIndex, mPositions.size() - 1);
    }

    private boolean inRegion(P current, int start, int end)
    {
        return current.Position >= start && current.Position <= end;
    }

    @NotNull
    private P current()
    {
        return mPositions.get(mIndex);
    }
}
