package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import jakarta.validation.constraints.NotNull;

public class CurationApplier
{
    @NotNull
    public static FindingRecord applyCurations(@NotNull FindingRecord findingRecord, @NotNull CurationRecord curationRecord)
    {
        Map<String, DriverCuration> curationByKey = curationRecord.driverCurations().stream()
                .collect(Collectors.toMap(DriverCuration::findingKey, c -> c));

        return FindingRecordBuilder.builder(findingRecord)
                .somaticSmallVariants(applyToDriverList(findingRecord.somaticSmallVariants(), curationByKey))
                .germlineSmallVariants(applyToDriverList(findingRecord.germlineSmallVariants(), curationByKey))
                .somaticGainDeletions(applyToDriverList(findingRecord.somaticGainDeletions(), curationByKey))
                .germlineGainDeletions(applyToDriverList(findingRecord.germlineGainDeletions(), curationByKey))
                .somaticDisruptions(applyToDriverList(findingRecord.somaticDisruptions(), curationByKey))
                .germlineDisruptions(applyToDriverList(findingRecord.germlineDisruptions(), curationByKey))
                .fusions(applyToDriverList(findingRecord.fusions(), curationByKey))
                .viruses(applyToDriverList(findingRecord.viruses(), curationByKey))
                .build();
    }

    @NotNull
    private static <T extends Driver> DriverFindingList<T> applyToDriverList(
            @NotNull DriverFindingList<T> driverFindingList, @NotNull Map<String, DriverCuration> curationByKey)
    {
        List<T> updatedFindings = driverFindingList.findings().stream()
                .map(finding -> applyCuration(finding, curationByKey))
                .toList();
        return new DriverFindingList<>(driverFindingList.status(), updatedFindings);
    }

    @NotNull
    @SuppressWarnings("unchecked")
    private static <T extends Driver> T applyCuration(@NotNull T driver, @NotNull Map<String, DriverCuration> curationByKey)
    {
        DriverCuration curation = curationByKey.get(driver.findingKey());
        if(curation == null)
        {
            return driver;
        }

        ReportedStatus newStatus = curation.reportableOverride() ? ReportedStatus.REPORTED : ReportedStatus.CANDIDATE;
        DriverFields updatedDriver = DriverFieldsBuilder.builder()
                .findingKey(driver.findingKey())
                .driverSource(driver.driverSource())
                .reportedStatus(newStatus)
                .driverInterpretation(driver.driverInterpretation())
                .driverLikelihood(driver.driverLikelihood())
                .build();

        if(driver instanceof SmallVariant v)
        {
            return (T) SmallVariantBuilder.builder(v).driver(updatedDriver).build();
        }
        else if(driver instanceof GainDeletion g)
        {
            return (T) GainDeletionBuilder.builder(g).driver(updatedDriver).build();
        }
        else if(driver instanceof Disruption d)
        {
            return (T) DisruptionBuilder.builder(d).driver(updatedDriver).build();
        }
        else if(driver instanceof Fusion f)
        {
            return (T) FusionBuilder.builder(f).driver(updatedDriver).build();
        }
        else if(driver instanceof Virus v)
        {
            return (T) VirusBuilder.builder(v).driver(updatedDriver).build();
        }
        throw new IllegalArgumentException("Unknown driver type: " + driver.getClass());
    }
}
