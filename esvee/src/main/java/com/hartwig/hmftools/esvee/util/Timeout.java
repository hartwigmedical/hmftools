package com.hartwig.hmftools.esvee.util;

import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.esvee.JunctionProcessingException;
import com.hartwig.hmftools.esvee.SVAConfig;

public class Timeout
{
    private final boolean mEnabled;
    private final long mDeadlineTimeNanos;
    private final Map<Object, Object> mContext = new LinkedHashMap<>();

    public Timeout(final SVAConfig config, final long timeoutNanos)
    {
        this(config.timeoutsEnabled(), timeoutNanos);
    }

    public Timeout(final boolean enabled, final long timeoutNanos)
    {
        mEnabled = enabled;
        mDeadlineTimeNanos = System.nanoTime() + timeoutNanos;
    }

    public void addContext(final Object object)
    {
        //mContext.put(object, object);
    }

    public void addContext(final String key, final Object object)
    {
        mContext.put(key, object);
    }

    public void checkTimeout()
    {
        if (mEnabled && System.nanoTime() > mDeadlineTimeNanos)
            throw TimeoutException.create(mContext);
    }

    public void checkTimeout(final String contextKeyToAddOnThrow, final Object object)
    {
        if (mEnabled && System.nanoTime() > mDeadlineTimeNanos)
        {
            addContext(contextKeyToAddOnThrow, object);
            throw TimeoutException.create(mContext);
        }
    }

    public static class TimeoutException extends JunctionProcessingException
    {
        @SuppressWarnings("unused") // For debugging
        private final Map<Object, Object> mContext;

        public TimeoutException()
        {
            super("Timeout expired!");
            mContext = new HashMap<>();
        }

        public TimeoutException(final String message, final Map<Object, Object> context)
        {
            super(message);
            mContext = context;
        }

        public static TimeoutException create(final Map<Object, Object> context)
        {
            if(context.isEmpty())
                return new TimeoutException();
            else
            {
                final String message = "Timeout expired!\nContext:\n" + context.entrySet().stream()
                        .map(entry -> entry.getKey() == entry.getValue()
                                ? entry.getValue().toString()
                                : entry.getKey().toString() + ": " + entry.getValue().toString())
                        .collect(Collectors.joining("\n"));
                return new TimeoutException(message, context);
            }
        }
    }
}
