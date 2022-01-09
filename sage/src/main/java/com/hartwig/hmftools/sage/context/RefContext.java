package com.hartwig.hmftools.sage.context;

import java.util.Collection;
import java.util.Map;

import javax.annotation.Nullable;

import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.read.ReadContext;

public class RefContext implements GenomePosition
{
    public final String Sample;
    public final String Chromosome;
    public final int MaxDepth;
    public final int Position;
    
    private final Map<String,AltContext> mAlts;

    private int mRawDepth;

    public RefContext(final String sample, final String chromosome, int position, int maxDepth)
    {
        Sample = sample;
        Chromosome = chromosome;
        Position = position;
        MaxDepth = maxDepth;
        mAlts = Maps.newHashMap();
    }

    public Collection<AltContext> alts()
    {
        return mAlts.values();
    }

    public boolean reachedLimit()
    {
        return mRawDepth >= MaxDepth;
    }

    public void refRead(boolean sufficientMapQuality)
    {
        if(sufficientMapQuality)
            mRawDepth++;
    }

    public void processAltRead(
            final String ref, final String alt, int baseQuality, boolean sufficientMapQuality,
            int numberOfEvents, final ReadContext readContext)
    {
        final AltContext altContext = getOrCreateAltContext(ref, alt);
        altContext.incrementAltRead(baseQuality);

        if(sufficientMapQuality)
            mRawDepth++;

        if(readContext != null && !readContext.incompleteCore())
        {
            altContext.addReadContext(numberOfEvents, readContext);
        }
    }

    @Override
    public String chromosome()
    {
        return Chromosome;
    }

    @Override
    public int position()
    {
        return Position;
    }

    public int rawDepth()
    {
        return mRawDepth;
    }

    @Override
    public boolean equals(@Nullable Object another)
    {
        if(this == another)
            return true;

        return another instanceof RefContext && equalTo((RefContext) another);
    }

    private boolean equalTo(RefContext another)
    {
        return Chromosome.equals(another.Chromosome) && Position == another.Position;
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Ints.hashCode(position());
        return h;
    }

    private AltContext getOrCreateAltContext(final String ref, final String alt)
    {
        final String refAltKey = ref + "|" + alt;
        AltContext altContext = mAlts.get(refAltKey);

        if(altContext == null)
        {
            altContext = new AltContext(this, ref, alt);
            mAlts.put(refAltKey, altContext);
        }

        return altContext;
    }
}
