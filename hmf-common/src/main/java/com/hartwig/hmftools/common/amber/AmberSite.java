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

    private boolean mSnpCheck;

    public AmberSite(final String chromosome, final int position, final String ref, final String alt, final boolean snpCheck)
    {
        Chromosome = chromosome;
        Position = position;
        Ref = ref;
        Alt = alt;
        mSnpCheck = snpCheck;
    }

    @Override
    public String chromosome() { return Chromosome; }

    @Override
    public int position() { return Position; }

    // convenience
    public String ref() { return Ref; }
    public String alt() { return Alt; }
    public boolean snpCheck() { return mSnpCheck; };
    public void setSnpCheck(boolean value) { mSnpCheck = value; };

    public String toString() { return format("%s:%d %s>%s %s", Chromosome, Position, Ref, Alt, mSnpCheck ? "snpcheck" : ""); }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof AmberSite && matches((AmberSite) another);
    }

    public boolean matches(final AmberSite another)
    {
        return matches(another.Chromosome, another.Position, another.Ref, another.Alt);
    }

    public boolean matches(final String chromosome, final int position, final String ref, final String alt)
    {
        return Chromosome.equals(chromosome) && Position == position && Ref.equals(ref) && Alt.equals(alt);
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + Chromosome.hashCode();
        h += (h << 5) + Position;
        h += (h << 5) + Ref.hashCode();
        h += (h << 5) + Alt.hashCode();
        h += (h << 5) + Booleans.hashCode(mSnpCheck);
        return h;
    }
}
