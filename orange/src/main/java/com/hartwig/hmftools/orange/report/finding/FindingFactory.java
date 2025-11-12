package com.hartwig.hmftools.orange.report.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.finding.FindingRecord;
import com.hartwig.hmftools.datamodel.finding.ImmutableFindingRecord;
import com.hartwig.hmftools.datamodel.finding.ImmutablePredictedTumorOrigin;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class FindingFactory {
    private FindingFactory() {
    }

    @NotNull
    public static FindingRecord create(@NotNull PurpleRecord purple,
            @NotNull LinxRecord linx,
            @Nullable VirusInterpreterData virusInterpreterData,
            @Nullable CuppaData cuppaData)
    {

        boolean containsTumorCells = purple.fit().qc().status().contains(PurpleQCStatus.FAIL_NO_TUMOR);

        ImmutableFindingRecord.Builder builder = ImmutableFindingRecord.builder();

        builder.allSomaticVariants(SmallVariantFactory.create(purple.allSomaticVariants(), purple.somaticDrivers()));

        List<PurpleVariant> allGermlineVariants = purple.allGermlineVariants();
        List<PurpleDriver> germlineDrivers = purple.germlineDrivers();
        if(allGermlineVariants != null && germlineDrivers != null) {
            builder.allGermlineVariants(SmallVariantFactory.create(allGermlineVariants,
                    germlineDrivers));
        }


        builder.allSomaticCopyNumbers(CopyNumberFactory.convert(purple.allSomaticGainsDels()));

        builder.allSomaticFusions(FusionFactory.convert(linx.allSomaticFusions()));

        builder.allSomaticDisruptions(DisruptionFactory.convert(linx.allSomaticBreakends(),
                linx.allSomaticStructuralVariants(),
                containsTumorCells));

        if(virusInterpreterData != null)
        {
            builder.allViruses(VirusFactory.convert(virusInterpreterData.allViruses()));
        }

        if(cuppaData != null)
        {
            builder.predictedTumorOrigin(ImmutablePredictedTumorOrigin.builder()
                            .findingKey(FindingKeys.findingKey(cuppaData.bestPrediction()))
                            .cancerType(cuppaData.bestPrediction().cancerType())
                            .likelihood(cuppaData.bestPrediction().likelihood())
                            .build());
        }

        //return report;
        return builder.build();
    }
}
