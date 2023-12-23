package com.hartwig.hmftools.svassembly.util;

import java.util.ArrayList;
import java.util.List;

public class FastFormat // FIXME: Do this and use for html page generation
{
    public static FastFormat compile(final String pattern)
    {
        final List<Appender> appenders = new ArrayList<>();

        final StringBuilder currentArg = new StringBuilder();
        for (int i = 0; i < pattern.length(); i++)
        {
            // The full format specifier looks like this:
            // %[argument_index$][flags][width][.precision][t]conversion
            // We only support %s for now :)


        }

        return null;
    }

    private static class RawArgAppender implements Appender
    {
        private final int mIndex;

        public RawArgAppender(final int index)
        {
            mIndex = index;
        }

        @Override
        public void append(final StringBuilder sb, final Object[] args)
        {
            sb.append(args[mIndex]);
        }
    }

    private static class ConstAppender implements Appender
    {
        private final String mValue;

        public ConstAppender(final String value)
        {
            mValue = value;
        }

        @Override
        public void append(final StringBuilder sb, final Object[] args)
        {
            sb.append(mValue);
        }
    }

    private interface Appender {
        void append(final StringBuilder sb, Object[] args);
    }
}
