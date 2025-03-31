package com.hartwig.hmftools.redux.common;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.sequencing.SequencingType;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

public class DuplicateGroupCollapseConfig
{
    public final SequencingType Sequencing;
    public final int SbxMaxDuplicateDistance;

    private static final String SBX_MAX_DUPLICATE_DISTANCE = "sbx_max_dup_dist";

    private static final String SBX_MAX_DUPLICATE_DISTANCE_DESC = "Max distance between the end of SBX fragments to be declared duplicates";

    public DuplicateGroupCollapseConfig(final SequencingType sequencing, final int sbxMaxDuplicateDistance)
    {
        Sequencing = sequencing;
        SbxMaxDuplicateDistance = sbxMaxDuplicateDistance;
    }

    @VisibleForTesting
    public DuplicateGroupCollapseConfig(final SequencingType sequencing)
    {
        this(sequencing, 0);
    }

    public static void addConfig(final ConfigBuilder configBuilder)
    {
        // sequencing is already defined in ReduxConfig

        configBuilder.addInteger(SBX_MAX_DUPLICATE_DISTANCE, SBX_MAX_DUPLICATE_DISTANCE_DESC, 0);
    }

    public static DuplicateGroupCollapseConfig from(final SequencingType sequencing, final ConfigBuilder configBuilder)
    {
        int sbxMaxDuplicateDistance = configBuilder.getInteger(SBX_MAX_DUPLICATE_DISTANCE);
        return new DuplicateGroupCollapseConfig(sequencing, sbxMaxDuplicateDistance);
    }
}
