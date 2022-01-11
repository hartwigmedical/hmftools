package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.SageVariant;

public abstract class BufferedPostProcessor implements Consumer<SageVariant>
{
    private final int mMaxDistance;

    private final ArrayDeque<SageVariant> mVariantQueue = new ArrayDeque<>();
    private final Consumer<SageVariant> mConsumer;

    public BufferedPostProcessor(int maxDistance, final Consumer<SageVariant> consumer)
    {
        mMaxDistance = maxDistance;
        mConsumer = consumer;
    }

    @Override
    public void accept(final SageVariant newVariant)
    {
        flush(newVariant);

        if(!mVariantQueue.isEmpty())
            processSageVariant(newVariant, mVariantQueue);
        
        mVariantQueue.add(newVariant);
    }

    protected abstract void processSageVariant(final SageVariant newVariant, final Collection<SageVariant> buffer);

    public void flush()
    {
        preFlush(mVariantQueue);
        mVariantQueue.forEach(mConsumer);
        mVariantQueue.clear();
    }

    protected void flush(final SageVariant variant)
    {
        if(mVariantQueue.isEmpty())
            return;

        final List<SageVariant> flushed = Lists.newArrayList();
        final Iterator<SageVariant> iterator = mVariantQueue.iterator();

        while(iterator.hasNext())
        {
            final SageVariant queuedVariant = iterator.next();

            int entryEnd = queuedVariant.position() + queuedVariant.ref().length() - 1;
            if(!queuedVariant.chromosome().equals(variant.chromosome()) || entryEnd < variant.position() - mMaxDistance)
            {
                iterator.remove();
                flushed.add(queuedVariant);
            }
            else
            {
                break;
            }
        }

        if(flushed.isEmpty())
            return;

        preFlush(flushed);
        flushed.forEach(mConsumer);
        flushed.clear();
    }

    protected void preFlush(final Collection<SageVariant> variants) { }

    public static boolean longerContainsShorter(final SageVariant shorter, final SageVariant longer)
    {
        return longerContainsShorter(shorter.variant(), longer.variant());
    }

    public static boolean longerContainsShorter(final VariantHotspot shorter, final VariantHotspot longer)
    {
        int longerStart = longer.position();
        int longerEnd = longer.end();

        int shorterStart = shorter.position();
        int shorterEnd = shorter.end();

        if(shorterStart < longerStart || shorterEnd > longerEnd)
            return false;

        final String shorterAlt = shorter.alt();

        int offset = shorterStart - longerStart;
        final String longerAlt = new String(longer.alt().getBytes(), offset, shorter.alt().length());
        return shorterAlt.equals(longerAlt);
    }

}
