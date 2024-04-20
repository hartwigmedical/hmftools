package com.hartwig.hmftools.esvee.prep.types;

import com.hartwig.hmftools.esvee.prep.PrepConfig;

public class JunctionsConfig
{
    public final ReadFilters ReadFiltering;
    public final int ReadLength;
    public final boolean TrimReadId;
    public final boolean UnpairedReads;
    public final boolean TrackRemotes;
    public final boolean PerfDebug;

    public JunctionsConfig(
            final ReadFilters readFiltering, final int readLength, final boolean trimReadId,
            final boolean unpairedReads, final boolean trackRemotes, final boolean perfDebug)
    {
        ReadFiltering = readFiltering;
        ReadLength = readLength;
        TrimReadId = trimReadId;
        UnpairedReads = unpairedReads;
        TrackRemotes = trackRemotes;
        PerfDebug = perfDebug;
    }

    public static JunctionsConfig from(final PrepConfig svConfig)
    {
        return new JunctionsConfig(
                svConfig.ReadFiltering, svConfig.ReadLength, svConfig.TrimReadId, svConfig.UnpairedReads,
                svConfig.TrackRemotes, svConfig.PerfDebug);
    }
}
