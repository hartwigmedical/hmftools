package com.hartwig.hmftools.orange.report.finding;

import java.util.Collection;
import java.util.List;

import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.Fusion;
import com.hartwig.hmftools.datamodel.finding.ImmutableFusion;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;

import org.jetbrains.annotations.NotNull;

public class FusionFactory {

    private FusionFactory() {
    }

    @NotNull
    public static List<Fusion> convert(@NotNull Collection<LinxFusion> linxFusions)
    {
        return linxFusions.stream().map(FusionFactory::convert)
                .toList();
    }

    private static Fusion convert(LinxFusion linxFusion)
    {
        return ImmutableFusion.builder().from(linxFusion)
                .findingKey(FindingKeys.findingKey(linxFusion))
                .isReportable(linxFusion.reported())
                .isCandidate(false)
                .driverInterpretation(convert(linxFusion.driverLikelihood()))
                .build();
    }

    private static DriverInterpretation convert(FusionLikelihoodType fusionLikelihoodType)
    {
        return switch (fusionLikelihoodType)
        {
            case HIGH -> DriverInterpretation.HIGH;
            case LOW, NA -> DriverInterpretation.LOW;
        };
    }
}
