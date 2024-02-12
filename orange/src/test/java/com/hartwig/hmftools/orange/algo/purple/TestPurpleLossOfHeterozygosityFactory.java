package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.purple.GeneProportion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleLossOfHeterozygosity;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleLossOfHeterozygosityFactory
{
    private TestPurpleLossOfHeterozygosityFactory()
    {
    }

    @NotNull
    public static ImmutablePurpleLossOfHeterozygosity.Builder builder()
    {
        return ImmutablePurpleLossOfHeterozygosity.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(true)
                .geneProportion(GeneProportion.PARTIAL_GENE)
                .minCopies(0)
                .maxCopies(0);
    }
}
