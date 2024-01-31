package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleGeneCopyNumberFactory
{
    private TestPurpleGeneCopyNumberFactory()
    {
    }

    @NotNull
    public static ImmutablePurpleGeneCopyNumber.Builder builder()
    {
        return ImmutablePurpleGeneCopyNumber.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .geneName(Strings.EMPTY)
                .minCopyNumber(0)
                .minMinorAlleleCopyNumber(0);
    }
}
