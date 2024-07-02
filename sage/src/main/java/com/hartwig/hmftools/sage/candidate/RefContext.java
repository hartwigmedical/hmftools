package com.hartwig.hmftools.sage.candidate;

import java.util.Collection;
import java.util.Map;

import com.google.common.collect.Maps;
import com.google.common.primitives.Ints;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.sage.common.RefSequence;
import com.hartwig.hmftools.sage.common.VariantReadContextBuilder;

import htsjdk.samtools.SAMRecord;

public class RefContext extends BasePosition
{
    private Map<String,AltContext> mAlts;

    public RefContext(final String chromosome, int position)
    {
        super(chromosome, position);
        mAlts = null;
    }

    public Collection<AltContext> altContexts()
    {
        return mAlts != null ? mAlts.values() : null;
    }

    public void processAltRead(
            final String ref, final String alt, int numberOfEvents, final SAMRecord read, final int variantReadIndex,
            final VariantReadContextBuilder readContextBuilder, final RefSequence refSequence)
    {
        final AltContext altContext = getOrCreateAltContext(ref, alt);
        altContext.incrementAltRead();

        altContext.addReadContext(numberOfEvents, read, variantReadIndex, readContextBuilder, refSequence);
    }

    public String chromosome()
    {
        return Chromosome;
    }
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

    private boolean equalTo(final RefContext another) { return matches(another); }

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
