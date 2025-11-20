package com.hartwig.hmftools.orange.algo.purple;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableGainDeletion;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.orange.report.finding.FindingKeys;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestGainDeletionFactory
{
    @NotNull
    public static GainDeletion createGainDel(@NotNull String gene, @NotNull CopyNumberInterpretation interpretation)
    {
        return builder()
                .findingKey(String.format("gainDel[%s]", gene))
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.HIGH)
                .gene(gene).interpretation(interpretation).build();
    }

    @NotNull
    public static ImmutableGainDeletion.Builder builder()
    {
        return ImmutableGainDeletion.builder()
                .findingKey(Strings.EMPTY)
                .reportedStatus(ReportedStatus.REPORTED)
                .driverInterpretation(DriverInterpretation.HIGH)
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
