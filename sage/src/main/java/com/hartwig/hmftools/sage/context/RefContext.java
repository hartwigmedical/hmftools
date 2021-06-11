package com.hartwig.hmftools.sage.context;

import java.util.Collection;
import java.util.HashMap;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.primitives.Longs;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.read.ReadContext;

import org.jetbrains.annotations.NotNull;

public class RefContext implements GenomePosition
{
    public final String Sample;
    public final String Chromosome;
    public final int MaxDepth;
    public final long Position;
    
    private final Map<String, AltContext> mAlts;

    private int rawDepth;

    public RefContext(final String sample, final String chromosome, final long position, final int maxDepth)
    {
        Sample = sample;
        Chromosome = chromosome;
        Position = position;
        MaxDepth = maxDepth;
        mAlts = new HashMap<>();
    }

    @NotNull
    public Collection<AltContext> alts()
    {
        return mAlts.values();
    }

    public boolean reachedLimit()
    {
        return rawDepth >= MaxDepth;
    }

    public void refRead(boolean sufficientMapQuality)
    {
        if(sufficientMapQuality)
        {
            rawDepth++;
        }
    }

    public void altRead(@NotNull final String ref, @NotNull final String alt, int baseQuality, boolean sufficientMapQuality,
            int numberOfEvents,
            @Nullable final ReadContext readContext)
    {
        final AltContext altContext = altContext(ref, alt);
        altContext.incrementAltRead(baseQuality);
        if(sufficientMapQuality)
        {
            rawDepth++;
        }

        if(readContext != null && !readContext.incompleteCore())
        {
            altContext.addReadContext(numberOfEvents, readContext);
        }
    }

    @NotNull
    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public long position()
    {
        return Position;
    }

    public int rawDepth()
    {
        return rawDepth;
    }

    @NotNull
    public String sample()
    {
        return Sample;
    }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
        {
            return true;
        }
        return another instanceof RefContext && equalTo((RefContext) another);
    }

    private boolean equalTo(RefContext another)
    {
        return chromosome().equals(another.chromosome()) && position() == another.position();
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Longs.hashCode(position());
        return h;
    }

    @NotNull
    private AltContext altContext(@NotNull final String ref, @NotNull final String alt)
    {
        final String refAltKey = ref + "|" + alt;
        return mAlts.computeIfAbsent(refAltKey, key -> new AltContext(this, ref, alt));
    }
}
