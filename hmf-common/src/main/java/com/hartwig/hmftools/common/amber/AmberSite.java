package com.hartwig.hmftools.common.amber;

import static java.lang.String.format;

import javax.annotation.Nullable;

import com.google.common.primitives.Booleans;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

public class AmberSite implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    public final String Ref;
    public final String Alt;
    public final boolean SnpCheck;

    public AmberSite(final String chromosome, final int position, final String ref, final String alt, final boolean snpCheck)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        SnpCheck = snpCheck;
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int position() { return Position; }

    // convenience
    public String ref() { return Ref; }
    public String alt() { return Alt; }
    public boolean snpCheck() { return SnpCheck; };

    public String toString() { return format("%s:%d %s>%s %s", Chromosome, Position, Ref, Alt, SnpCheck ? "snpcheck" : ""); }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof AmberSite && equalTo((AmberSite) another);
    }

    private boolean equalTo(final AmberSite another)
    {
        return Chromosome.equals(another.Chromosome)
                && Position == another.Position
                && Ref.equals(another.Ref)
                && Alt.equals(another.Alt)
                && SnpCheck == another.SnpCheck;
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + Chromosome.hashCode();
        h += (h << 5) + Position;
        h += (h << 5) + Ref.hashCode();
        h += (h << 5) + Alt.hashCode();
        h += (h << 5) + Booleans.hashCode(SnpCheck);
        return h;
    }
}
