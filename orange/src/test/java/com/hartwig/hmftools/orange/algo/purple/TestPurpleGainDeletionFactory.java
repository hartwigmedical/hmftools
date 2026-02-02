package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleGainDeletionFactory
{
    @NotNull
    public static PurpleGainDeletion createGainDel(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation)
    {
        return builder(gene).interpretation(interpretation).build();
    }

    @NotNull
    public static ImmutablePurpleGainDeletion.Builder builder()
    {
        return ImmutablePurpleGainDeletion.builder()
                .driver(driverBuilder().build())
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .minCopies(0)
                .maxCopies(0)
                .minMinorAlleleCopies(0);
    }

    @NotNull
    public static ImmutablePurpleGainDeletion.Builder builder(@NotNull String gene)
    {
        return builder().driver(driverBuilder().gene(gene).build());
    }

    @NotNull
    public static ImmutablePurpleDriver.Builder driverBuilder()
    {
        return ImmutablePurpleDriver.builder()
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .type(PurpleDriverType.DEL)
                .driverLikelihood(0.0)
                .likelihoodMethod(PurpleLikelihoodMethod.DEL)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.HIGH)
                .isCanonical(true);
    }
}
