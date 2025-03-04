package com.hartwig.hmftools.pave.reverse;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public final class ChangeLocation
{
    @NotNull
    public final String mChromosome;

    public final int Location;

    public ChangeLocation(@NotNull final String chromosome, final int location)
    {
        this.mChromosome = chromosome;
        Location = location;
    }
    @Override
    public String toString()
    {
        return "ChangeLocation{" +
                "chromosome=" + mChromosome +
                ", location=" + Location +
                '}';
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final ChangeLocation that = (ChangeLocation) o;
        return Location == that.Location && Objects.equals(mChromosome, that.mChromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Location);
    }
}
