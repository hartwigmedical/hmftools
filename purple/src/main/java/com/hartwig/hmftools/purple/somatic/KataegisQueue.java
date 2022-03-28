package com.hartwig.hmftools.purple.somatic;

import static com.hartwig.hmftools.common.variant.SomaticVariantFactory.KATAEGIS_FLAG;
import static com.hartwig.hmftools.purple.PurpleCommon.PPL_LOGGER;

import java.util.ArrayDeque;
import java.util.Deque;
import java.util.function.Consumer;
import java.util.function.Predicate;

public class KataegisQueue
{
    static final long MAX_ABS_DISTANCE = 2000;
    static final long MAX_AVG_DISTANCE = 1000;
    static final int MIN_COUNT = 3;

    private final Predicate<SomaticVariant> mCandidate;

    private final Consumer<SomaticVariant> mConsumer;
    private final Deque<SomaticVariant> mBuffer;
    private final String mIdPrefix;

    private int mIdentifier;

    public KataegisQueue(final String idPrefix, final Predicate<SomaticVariant> candidate, final Consumer<SomaticVariant> consumer)
    {
        mCandidate = candidate;
        mConsumer = consumer;
        mIdPrefix = idPrefix;
        mBuffer = new ArrayDeque<>();
        mIdentifier = 0;
    }

    public void processVariant(final SomaticVariant context)
    {
        if(!mBuffer.isEmpty())
        {
            final SomaticVariant previous = mBuffer.peekLast();
            if(context.position() - previous.position() > MAX_ABS_DISTANCE || !previous.chromosome().equals(context.chromosome()))
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
            final SomaticVariant first = mBuffer.peekFirst();

            if(!mCandidate.test(first))
            {
                SomaticVariant var = mBuffer.pollFirst();
                if(mConsumer != null)
                    mConsumer.accept(var);
            }
            else
            {
                final KataegisWindow window = longestViableWindow(first);
                final boolean isWindowViable = window.isViable(MIN_COUNT, MAX_AVG_DISTANCE);
                if(isWindowViable)
                {
                    mIdentifier++;
                }

                while(!mBuffer.isEmpty())
                {
                    final SomaticVariant peek = mBuffer.peekFirst();
                    if(peek.position() > window.end())
                    {
                        return;
                    }

                    if(isWindowViable && mCandidate.test(peek))
                    {
                        peek.context().getCommonInfo().putAttribute(KATAEGIS_FLAG, mIdPrefix + "_" + mIdentifier, true);
                        PPL_LOGGER.debug("var({}) added kataegis flag: {}", peek, peek.context().getAttribute(KATAEGIS_FLAG));
                    }

                    SomaticVariant var = mBuffer.pollFirst();
                    if(mConsumer != null)
                        mConsumer.accept(var);
                }
            }
        }
    }

    private KataegisWindow longestViableWindow(final SomaticVariant first)
    {
        KataegisWindow result = new KataegisWindow(first);
        final KataegisWindow window = new KataegisWindow(first);

        for(SomaticVariant variant : mBuffer)
        {
            if(variant.position() - window.end() > MAX_ABS_DISTANCE)
            {
                return result;
            }

            if(mCandidate.test(variant))
            {
                window.add(variant);
            }

            if(window.isViable(MIN_COUNT, MAX_AVG_DISTANCE))
            {
                result = new KataegisWindow(window);
            }
        }

        return result;
    }
}
