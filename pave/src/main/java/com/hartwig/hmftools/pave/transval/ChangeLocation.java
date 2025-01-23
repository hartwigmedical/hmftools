package com.hartwig.hmftools.pave.transval;

import java.util.Objects;

import org.jetbrains.annotations.NotNull;

public final class ChangeLocation
{
    @NotNull
    public final String Chromosome;

    public final int Location;

    public ChangeLocation(@NotNull final String chromosome, final int location)
    {
        this.Chromosome = chromosome;
        Location = location;
    }
    @Override
    public String toString()
    {
        return "ChangeLocation{" +
                "chromosome=" + Chromosome +
                ", Location=" + Location +
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
        return Location == that.Location && Objects.equals(Chromosome, that.Chromosome);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(Chromosome, Location);
    }
}
