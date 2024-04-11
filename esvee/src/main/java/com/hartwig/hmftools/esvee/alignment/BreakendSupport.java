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

    public String toString()
    {
        return format("support(split=%d disc=%d) strand(fwd=%d rev=%d)",
                SplitFragments, DiscordantFragments, ForwardReads, ReverseReads);
    }


    public static List<BreakendSupport> buildSampleData(final JunctionAssembly assembly, final List<String> sampleNames)
    {
        // count final support for each breakend by now testing each support read vs the alignment position
        // for reads from other locations, test vs the junction position in the assembly
        List<BreakendSupport> sampleDataList = Lists.newArrayListWithExpectedSize(sampleNames.size());

        sampleNames.forEach(x -> sampleDataList.add(new BreakendSupport()));

        int currentSampleIndex = -1;
        BreakendSupport sampleData = null;

        for(SupportRead support : assembly.support())
        {
            if(support.type() == SupportType.JUNCTION_MATE)
                continue;

            int sampleIndex = support.sampleIndex();
            String readSampleId = sampleNames.get(sampleIndex);

            if(currentSampleIndex != sampleIndex)
            {
                currentSampleIndex = sampleIndex;
                sampleData = sampleDataList.get(currentSampleIndex);
            }

            if(support.type().isSplitSupport())
                ++sampleData.SplitFragments;
            else
                ++sampleData.DiscordantFragments;
        }

        return sampleDataList;
    }

    /*
    private int getSampleIdIndex(final String sampleId)
    {
        for(int i = 0; i < mSampleNames.size(); ++i)
        {
            if(mSampleNames.get(i).equals(sampleId))
                return i;
        }

        return -1;
    }
    */
}
