package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestPurpleGainDeletionFactory
{
    public static PurpleGainDeletion createGainDel(final String gene)
    {
        return builder(gene).build();
    }

    public static ImmutablePurpleGainDeletion.Builder builder()
    {
        return ImmutablePurpleGainDeletion.builder()
                .driver(driverBuilder().build())
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .minCopyNumber(0)
                .maxCopyNumber(0)
                .relativeCopyNumber(0)
                .exonRange("")
                .minMinorAlleleCopies(0)
                .tpm(null)
                .tpmPercentile(null)
                .tpmFoldChange(null);
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
