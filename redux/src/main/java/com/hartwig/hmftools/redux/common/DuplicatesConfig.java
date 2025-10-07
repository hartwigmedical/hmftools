package com.hartwig.hmftools.redux.common;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class DuplicatesConfig
{
    public final int SbxMaxDuplicateDistance;

    private static final String SBX_MAX_DUPLICATE_DISTANCE = "sbx_max_dup_dist";
    private static final String SBX_MAX_DUPLICATE_DISTANCE_DESC = "Max distance between the end of SBX fragments to be declared duplicates";

    public DuplicatesConfig(final int sbxMaxDuplicateDistance)
    {
        SbxMaxDuplicateDistance = sbxMaxDuplicateDistance;
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        // sequencing is already defined in ReduxConfig
        configBuilder.addInteger(SBX_MAX_DUPLICATE_DISTANCE, SBX_MAX_DUPLICATE_DISTANCE_DESC, 0);
    }

    public static DuplicatesConfig from(final ConfigBuilder configBuilder)
    {
        int sbxMaxDuplicateDistance = configBuilder.getInteger(SBX_MAX_DUPLICATE_DISTANCE);
        return new DuplicatesConfig(sbxMaxDuplicateDistance);
    }
}
