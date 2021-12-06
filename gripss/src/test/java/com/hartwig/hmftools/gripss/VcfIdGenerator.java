package com.hartwig.hmftools.gripss;

public class VcfIdGenerator
{
    private int mEventId;

    private static final char SUFFIX_FIRST = 'o';
    private static final char SUFFIX_SECOND = 'h';

    public VcfIdGenerator()
    {
        mEventId = 0;
    }

    public String nextVcfId(boolean isFirst)
    {
        ++mEventId;
        return vcfId(suffix(isFirst));
    }

    public String nextEventId()
    {
        ++mEventId;
        return eventId();
    }

    public static String vcfId(final String eventId, boolean isFirst)
    {
        return eventId + (isFirst ? SUFFIX_FIRST : SUFFIX_SECOND);
    }

    // public String current(boolean isFirst) { return formatId(suffix(isFirst)); }

    private String eventId() { return String.format("vid_%d", mEventId); }

    private String vcfId(final char suffix) { return eventId() + suffix; }

    private char suffix(boolean isFirst) { return isFirst ? SUFFIX_FIRST : SUFFIX_SECOND; }
}
