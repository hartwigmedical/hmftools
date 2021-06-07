package com.hartwig.hmftools.sage.phase;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.function.Consumer;

import com.hartwig.hmftools.common.genome.position.GenomePosition;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.variant.SageVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public abstract class BufferedPostProcessor implements Consumer<SageVariant>
{

    private final int maxDistance;

    private final ArrayDeque<SageVariant> buffer = new ArrayDeque<>();
    private final Consumer<SageVariant> consumer;

    public BufferedPostProcessor(int maxDistance, final Consumer<SageVariant> consumer)
    {
        this.maxDistance = maxDistance;
        this.consumer = consumer;
    }

    @Override
    public void accept(final SageVariant newVariant)
    {
        flush(newVariant);
        processSageVariant(newVariant, buffer);
        buffer.add(newVariant);
    }

    protected abstract void processSageVariant(@NotNull final SageVariant newVariant, @NotNull final Collection<SageVariant> buffer);

    public static boolean longerContainsShorter(@NotNull final SageVariant shorter, @NotNull final SageVariant longer)
    {
        return longerContainsShorter(shorter.variant(), longer.variant());
    }

    public static boolean longerContainsShorter(@NotNull final VariantHotspot shorter, @NotNull final VariantHotspot longer)
    {
        long longerStart = longer.position();
        long longerEnd = longer.end();

        long shorterStart = shorter.position();
        long shorterEnd = shorter.end();

        if(shorterStart < longerStart || shorterEnd > longerEnd)
        {
            return false;
        }

        final String shorterAlt = shorter.alt();

        int offset = (int) (shorterStart - longerStart);
        final String longerAlt = new String(longer.alt().getBytes(), offset, shorter.alt().length());
        return shorterAlt.equals(longerAlt);
    }

    public void flush()
    {
        preFlush(buffer);
        buffer.forEach(consumer);
        buffer.clear();
    }

    protected void flush(@NotNull final GenomePosition position)
    {
        final List<SageVariant> flushed = Lists.newArrayList();
        final Iterator<SageVariant> iterator = buffer.iterator();
        while(iterator.hasNext())
        {
            final SageVariant entry = iterator.next();
            long entryEnd = entry.position() + entry.ref().length() - 1;
            if(!entry.chromosome().equals(position.chromosome()) || entryEnd < position.position() - maxDistance)
            {
                iterator.remove();
                flushed.add(entry);
            }
            else
            {
                break;
            }
        }

        if(!flushed.isEmpty())
        {
            preFlush(flushed);
        }

        flushed.forEach(consumer);
        flushed.clear();
    }

    protected void preFlush(@NotNull final Collection<SageVariant> variants)
    {

    }
}
