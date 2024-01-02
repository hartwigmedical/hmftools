package com.hartwig.hmftools.sage.candidate_;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.genome.position.GenomePosition;

import org.jetbrains.annotations.Nullable;

/**
 * A collection of AltContexts at a GenomePosition that we use to process alt read contexts at this GenomePosition.
 */
public class RefContext_ implements GenomePosition
{
    public final String Chromosome;
    public final int Position;

    private Map<String, AltContext_> mAlts;

    public RefContext_(final String chromosome, int position)
    {
        Chromosome = chromosome;
        Position = position;
        mAlts = null;
    }

    private AltContext_ getOrCreateAltContext(final String ref, final String alt)
    {
        if (mAlts == null)
            mAlts = Maps.newHashMap();

        String key = ref + "|" + alt;
        AltContext_ altContext = mAlts.get(key);
        if (altContext == null)
        {
            altContext = new AltContext_(this, ref, alt);
            mAlts.put(key, altContext);
        }

        return altContext;
    }

    /**
     * Gets the alt context, register the read with the baseQuality, and then adds the read context to the alt context using numberOfEvents.
     */
    public void processAltRead(final String ref, final String alt, int baseQuality, int numberOfEvents, @Nullable final ReadContext_ readContext)
    {
        final AltContext_ altContext = getOrCreateAltContext(ref, alt);
        altContext.incrementAltRead(baseQuality);

        if(readContext != null && !readContext.hasIncompleteCore())
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

    public Collection<AltContext_> altContexts()
    {
        return mAlts != null ? mAlts.values() : null;
    }

    private boolean equalTo(RefContext_ another)
    {
        return Chromosome.equals(another.Chromosome) && Position == another.Position;
    }

    /**
     * This is only based on Chromosome and Position.
     */
    @Override
    public boolean equals(final Object another)
    {
        if(this == another)
            return true;

        return another instanceof RefContext_ && equalTo((RefContext_) another);
    }

    @Override
    public int hashCode()
    {
        int h = 5381;
        h += (h << 5) + chromosome().hashCode();
        h += (h << 5) + Ints.hashCode(position());
        return h;
    }
}
