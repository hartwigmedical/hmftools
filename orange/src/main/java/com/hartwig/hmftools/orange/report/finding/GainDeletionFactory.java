package com.hartwig.hmftools.orange.report.finding;

import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.datamodel.finding.GainDeletion;
import com.hartwig.hmftools.datamodel.finding.ImmutableCopyNumber;
import com.hartwig.hmftools.datamodel.purple.PurpleGainDeletion;

import org.jetbrains.annotations.NotNull;

public class GainDeletionFactory
{
    private GainDeletionFactory() {
    }

    @NotNull
    public static List<GainDeletion> convert(@NotNull Collection<DriverCatalog> purpleGainsDeletions)
    {
        return purpleGainsDeletions.stream().map(o -> (GainDeletion) ImmutableCopyNumber.builder()
                        .findingKey(FindingKeys.findingKey(o))
                        .reportedStatus(o.reportedStatus())
                        .driverInterpretation(o.driverInterpretation())
                        .build())
                .toList();
    }
}
