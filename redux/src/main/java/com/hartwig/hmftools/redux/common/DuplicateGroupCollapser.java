package com.hartwig.hmftools.redux.common;

import static com.hartwig.hmftools.common.sequencing.SequencingType.ULTIMA;

import java.util.List;

import com.hartwig.hmftools.common.sequencing.SequencingType;

import org.jetbrains.annotations.Nullable;

public abstract class DuplicateGroupCollapser
{
    public abstract FragmentCoordReads collapse(@Nullable final List<DuplicateGroup> duplicateGroups,
            @Nullable final List<ReadInfo> singleReads);

    public static DuplicateGroupCollapser fromSequencingType(final SequencingType sequencingType)
    {
        if(sequencingType == ULTIMA)
            return new UltimaDuplicateGroupCollapser();

        return null;
    }
}
