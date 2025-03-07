package com.hartwig.hmftools.pavereverse.serve;

import java.util.Objects;

import com.hartwig.hmftools.pavereverse.BaseSequenceChange;

import org.jetbrains.annotations.NotNull;

class FloatingHotspot
{
    @NotNull
    public final String mChromosome;

    @NotNull
    public final String Ref;

    @NotNull
    public final String Alt;

    public FloatingHotspot(BaseSequenceChange hotspot)
    {
        mChromosome = hotspot.mChromosome;
        Ref = hotspot.Ref;
        Alt = hotspot.Alt;
    }

    @Override
    public boolean equals(final Object o)
    {
        if(o == null || getClass() != o.getClass())
        {
            return false;
        }
        final FloatingHotspot that = (FloatingHotspot) o;
        return Objects.equals(mChromosome, that.mChromosome) && Objects.equals(Ref, that.Ref)
                && Objects.equals(Alt, that.Alt);
    }

    @Override
    public int hashCode()
    {
        return Objects.hash(mChromosome, Ref, Alt);
    }
}
