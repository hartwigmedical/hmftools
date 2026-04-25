package com.hartwig.hmftools.geneutils.panelfinder;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.ITEM_DELIM;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

public class HighDepthData extends ChrBaseRegion
{
    public final int SampleCount;
    public final int MinDepth;
    public final int MaxDepth;

    public HighDepthData(final ChrBaseRegion region, final int sampleCount, final int minDepth, final int maxDepth)
    {
        super(region.Chromosome, region.start(), region.end());

        SampleCount = sampleCount;
        MinDepth = minDepth;
        MaxDepth = maxDepth;
    }

    public String toString() { return format("%s: samples(%d) depth(%d-%d)", super.toString(), SampleCount, MinDepth, MaxDepth); }

    public static String toString(final List<HighDepthData> highDepths)
    {
        return highDepths.stream().map(x -> x.toString()).collect(Collectors.joining(ITEM_DELIM));
    }
}
