package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.Driver;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.purple.PurpleFit;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface FindingRecord
{
    @NotNull
    List<SmallVariant> driverSomaticSmallVariants();

    @NotNull
    default List<SmallVariant> driverSomaticSmallVariants(ReportedStatus reportedStatus)
    {
        return filterReported(driverSomaticSmallVariants(), reportedStatus);
    }

    @Nullable
    List<SmallVariant> driverGermlineSmallVariants();

    @Nullable
    default List<SmallVariant> driverGermlineSmallVariants(ReportedStatus reportedStatus)
    {
        return filterReported(driverGermlineSmallVariants(), reportedStatus);
    }

    @NotNull
    List<GainDeletion> driverSomaticGainDeletions();

    @NotNull
    default List<GainDeletion> driverSomaticGainDeletions(ReportedStatus reportedStatus)
    {
        return filterReported(driverSomaticGainDeletions(), reportedStatus);
    }

    @Nullable
    List<GainDeletion> driverGermlineGainDeletions();

    @NotNull
    List<Fusion> driverSomaticFusions();

    @NotNull
    default List<Fusion> driverSomaticFusions(ReportedStatus reportedStatus)
    {
        return filterReported(driverSomaticFusions(), reportedStatus);
    }

    @NotNull
    List<Disruption> driverSomaticDisruptions();

    @NotNull
    default List<Disruption> driverSomaticDisruptions(ReportedStatus reportedStatus)
    {
        return filterReported(driverSomaticDisruptions(), reportedStatus);
    }

    @Nullable
    List<Disruption> driverGermlineDisruptions();

    @Nullable
    default List<Disruption> driverGermlineDisruptions(ReportedStatus reportedStatus)
    {
        return filterReported(driverGermlineDisruptions(), reportedStatus);
    }

    @NotNull
    List<ChromosomeArmCopyNumber> chromosomeArmCopyNumbers();

    @Nullable
    List<Virus> driverViruses();

    @Nullable
    default List<Virus> driverViruses(ReportedStatus reportedStatus)
    {
        return filterReported(driverViruses(), reportedStatus);
    }

    @NotNull
    Findings<HlaAllele> hlaFindings();

    @Nullable
    MicrosatelliteStability microsatelliteStability();

    @Nullable
    TumorMutationStatus tumorMutationStatus();

    @Nullable
    PredictedTumorOrigin predictedTumorOrigin();

    @Nullable
    HomologousRecombination homologousRecombination();

    @NotNull
    OrangeRefGenomeVersion refGenomeVersion();

    @NotNull
    ExperimentType experimentType();

    @Nullable
    String pipelineVersion();

    @NotNull
    PurpleFit purpleFit();

    private static <T extends Driver> List<T> filterReported(List<T> drivers, ReportedStatus reportedStatus)
    {
        if (drivers == null) { return null; }
        return drivers.stream().filter(o -> o.reportedStatus() == reportedStatus).toList();
    }
}
