package com.hartwig.hmftools.bamtools.copyfastqtags;

import java.util.List;

public class FastqReadTags
{
    public final String ReadName;
    public final List<FastqTag> Tags;

    public FastqReadTags(final String readName, final List<FastqTag> tags)
    {
        ReadName = readName;
        Tags = tags;
    }
}
