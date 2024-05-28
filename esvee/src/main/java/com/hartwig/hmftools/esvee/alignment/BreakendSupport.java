package com.hartwig.hmftools.esvee.alignment;

import static java.lang.String.format;

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

    public double strandBias()
    {
        int totalReads = ForwardReads + ReverseReads;
        return totalReads > 0 ? ForwardReads / (double)totalReads : 0;
    }

    public String toString()
    {
        return format("support(split=%d disc=%d) strand(fwd=%d rev=%d)",
                SplitFragments, DiscordantFragments, ForwardReads, ReverseReads);
    }

}
