package com.hartwig.hmftools.orange.report.finding;

import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.datamodel.finding.CopyNumber;
import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ImmutableCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

import org.jetbrains.annotations.NotNull;

public class CopyNumberFactory
{
    private CopyNumberFactory() {
    }

    @NotNull
    public static List<CopyNumber> convert(@NotNull Collection<PurpleGainDeletion> purpleGainsDeletions)
    {
        return purpleGainsDeletions.stream().map(o -> (CopyNumber) ImmutableCopyNumber.builder().from(o)
                        .isReportable(true) // TODOHWL: fix this
                        .driverInterpretation(DriverInterpretation.HIGH) // TODOHWL: fix this
                        .build())
                .toList();
    }
}
