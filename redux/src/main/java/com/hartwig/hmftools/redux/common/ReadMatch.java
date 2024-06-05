package com.hartwig.hmftools.redux.common;

public class ReadMatch
{
    public final boolean Matched;
    public final FragmentStatus Status;

    public static final ReadMatch NO_READ_MATCH = new ReadMatch(false, null);

    public ReadMatch(final boolean matched, final FragmentStatus status)
    {
        Matched = matched;
        Status = status;
    }
}
