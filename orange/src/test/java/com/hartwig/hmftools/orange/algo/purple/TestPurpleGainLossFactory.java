package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainLoss;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDel;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleGainLossFactory
{
    @NotNull
    public static PurpleGainDel createGainLoss(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation)
    {
        return builder().gene(gene).interpretation(interpretation).build();
    }

    @NotNull
    public static ImmutablePurpleGainLoss.Builder builder()
    {
        return ImmutablePurpleGainLoss.builder()
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
