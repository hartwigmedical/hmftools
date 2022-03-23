package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.KATAEGIS_FLAG;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.function.Consumer;
import java.util.function.Predicate;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.VariantContext;

public class KataegisQueue
{
    static final long MAX_ABS_DISTANCE = 2000;
    static final long MAX_AVG_DISTANCE = 1000;
    static final int MIN_COUNT = 3;

    private final Predicate<VariantContext> mCandidate;

    @Nullable
    private final Consumer<VariantContext> mConsumer;
    private final Deque<VariantContext> mBuffer;
    private final String mIdPrefix;

    private int identifier = 0;

    public KataegisQueue(final String idPrefix, final Predicate<VariantContext> candidate, final Consumer<VariantContext> consumer)
    {
        mCandidate = candidate;
        mConsumer = consumer;
        mIdPrefix = idPrefix;
        mBuffer = new ArrayDeque<>();
    }

    public void accept(final VariantContext context)
    {
        if(!mBuffer.isEmpty())
        {
            final VariantContext previous = mBuffer.peekLast();
            if(context.getStart() - previous.getStart() > MAX_ABS_DISTANCE || !previous.getContig().equals(context.getContig()))
            {
                flush();
            }
        }

        mBuffer.add(context);
    }

    public void flush()
    {
        while(!mBuffer.isEmpty())
        {
            processFirstContext();
        }
    }

    private void processFirstContext()
    {
        if(!mBuffer.isEmpty())
        {
            final VariantContext first = mBuffer.peekFirst();

            if(!mCandidate.test(first))
            {
                if(mConsumer != null)
                    mConsumer.accept(mBuffer.pollFirst());
            }
            else
            {
                final KataegisWindow window = longestViableWindow(first);
                final boolean isWindowViable = window.isViable(MIN_COUNT, MAX_AVG_DISTANCE);
                if(isWindowViable)
                {
                    identifier++;
                }

                while(!mBuffer.isEmpty())
                {
                    final VariantContext peek = mBuffer.peekFirst();
                    if(peek.getStart() > window.end())
                    {
                        return;
                    }

                    if(isWindowViable && mCandidate.test(peek))
                    {
                        peek.getCommonInfo().putAttribute(KATAEGIS_FLAG, mIdPrefix + "_" + identifier, true);
                    }

                    if(mConsumer != null)
                        mConsumer.accept(mBuffer.pollFirst());
                }
            }
        }
    }

    private KataegisWindow longestViableWindow(final VariantContext first)
    {
        KataegisWindow result = new KataegisWindow(first);
        final KataegisWindow window = new KataegisWindow(first);

        for(VariantContext context : mBuffer)
        {
            if(context.getStart() - window.end() > MAX_ABS_DISTANCE)
            {
                return result;
            }

            if(mCandidate.test(context))
            {
                window.add(context);
            }

            if(window.isViable(MIN_COUNT, MAX_AVG_DISTANCE))
            {
                result = new KataegisWindow(window);
            }
        }

        return result;
    }
}
