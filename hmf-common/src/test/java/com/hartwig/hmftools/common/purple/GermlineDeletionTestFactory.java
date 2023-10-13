package com.hartwig.hmftools.common.purple;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class GermlineDeletionTestFactory
{
    @NotNull
    public static GermlineDeletion create(@NotNull String geneName)
    {
        return create(geneName, false, GermlineStatus.HET_DELETION, 0D, 0, 0);
    }

    @NotNull
    public static GermlineDeletion create(@NotNull String geneName, boolean reported, @NotNull GermlineStatus tumorStatus)
    {
        return create(geneName, reported, tumorStatus, 0D, 0, 0);
    }

    @NotNull
    public static GermlineDeletion create(@NotNull String geneName, boolean reported, @NotNull GermlineStatus tumorStatus,
            double tumorCopyNumber)
    {
        return create(geneName, reported, tumorStatus, tumorCopyNumber, 0, 0);
    }

    @NotNull
    public static GermlineDeletion create(@NotNull String geneName, boolean reported, @NotNull GermlineStatus tumorStatus,
            double tumorCopyNumber, int regionStart, int regionEnd)
    {
        return new GermlineDeletion(geneName,
                Strings.EMPTY,
                Strings.EMPTY,
                regionStart,
                regionEnd,
                0,
                0,
                0,
                GermlineDetectionMethod.SEGMENT,
                tumorStatus,
                tumorStatus,
                0D,
                tumorCopyNumber,
                Strings.EMPTY,
                0,
                reported);
    }
}
