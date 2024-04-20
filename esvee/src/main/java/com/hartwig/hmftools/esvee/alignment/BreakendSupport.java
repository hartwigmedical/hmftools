package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.esvee.assembly.output.VcfWriter;
import com.hartwig.hmftools.esvee.assembly.types.JunctionAssembly;
import com.hartwig.hmftools.esvee.assembly.types.SupportRead;
import com.hartwig.hmftools.esvee.assembly.types.SupportType;

public class BreakendSupport
{
    public int SplitFragments;
    public int DiscordantFragments;
    public int ForwardReads;
    public int ReverseReads;

    public BreakendSupport()
    {
        SplitFragments = 0;
        DiscordantFragments = 0;
        ForwardReads = 0;
        ReverseReads = 0;
    }

    public int totalSupport() { return SplitFragments + DiscordantFragments; }

    public String toString()
    {
        return format("support(split=%d disc=%d) strand(fwd=%d rev=%d)",
                SplitFragments, DiscordantFragments, ForwardReads, ReverseReads);
    }
}
