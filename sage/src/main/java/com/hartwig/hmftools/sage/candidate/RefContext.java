package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.sage.common.ReadContext;

public class RefContext implements GenomePosition
{
    public final String Chromosome;
    public final int Position;
    
    private Map<String,AltContext> mAlts;

    public RefContext(final String chromosome, int position)
    {
        Chromosome = chromosome;
        Position = position;
        mAlts = null;
    }

    public Collection<AltContext> altContexts()
    {
        return mAlts != null ? mAlts.values() : null;
    }

    public void processAltRead(final String ref, final String alt, int baseQuality, int numberOfEvents, final ReadContext readContext)
    {
        final AltContext altContext = getOrCreateAltContext(ref, alt);
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

    @Override
    public boolean equals(final Object another)
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
        if(mAlts == null)
            mAlts = Maps.newHashMap();

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
