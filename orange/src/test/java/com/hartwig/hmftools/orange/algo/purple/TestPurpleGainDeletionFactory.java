package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleGainDeletionFactory
{
    @NotNull
    public static PurpleGainDeletion createGainDel(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation)
    {
        return builder().gene(gene).interpretation(interpretation).build();
    }

    @NotNull
    public static ImmutablePurpleGainDeletion.Builder builder()
    {
        return ImmutablePurpleGainDeletion.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(0)
                .maxCopies(0);
    }
}
