package com.hartwig.hmftools.sage.select;

import java.util.List;
import java.util.Optional;
import java.util.function.Consumer;

import javax.annotation.concurrent.NotThreadSafe;

import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.NotNull;

@NotThreadSafe
public class PositionSelector<P extends GenomePosition>
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

        int currentCompare = Long.compare(current().position(), position);
        while(currentCompare >= 0 && mIndex > 0)
        {
            mIndex--;
            currentCompare = Long.compare(current().position(), position);
        }

        while(currentCompare < 0 && mIndex < mPositions.size() - 1)
        {
            mIndex++;
            currentCompare = Long.compare(current().position(), position);
        }

        return currentCompare == 0 ? Optional.of(current()) : Optional.empty();
    }

    public void select(final long start, final long end, final Consumer<P> handler)
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

    private boolean inRegion(P current, long start, long end)
    {
        return current.position() >= start && current.position() <= end;
    }

    @NotNull
    private P current()
    {
        return mPositions.get(mIndex);
    }
}
