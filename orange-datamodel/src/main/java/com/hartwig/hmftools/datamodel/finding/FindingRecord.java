package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

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
    List<SmallVariant> allSomaticVariants();

    @Gson.Ignore
    @NotNull
    default List<SmallVariant> reportableSomaticVariants()
    {
        return filterReported(allSomaticVariants());
    }

    @Nullable
    List<SmallVariant> allGermlineVariants();

    @Gson.Ignore
    @Nullable
    default List<SmallVariant> reportableGermlineVariants()
    {
        return filterReported(allGermlineVariants());
    }

    @NotNull
    List<CopyNumber> allSomaticCopyNumbers();

    @Gson.Ignore
    @NotNull
    default List<CopyNumber> reportableSomaticCopyNumbers()
    {
        return filterReported(allSomaticCopyNumbers());
    }

    @Nullable
    List<Fusion> allSomaticFusions();

    @Gson.Ignore
    @Nullable
    default List<Fusion> reportableSomaticFusions()
    {
        return filterReported(allSomaticFusions());
    }

    @Nullable
    List<Disruption> allSomaticDisruptions();

    @Gson.Ignore
    @Nullable
    default List<Disruption> reportableSomaticDisruptions()
    {
        return filterReported(allSomaticDisruptions());
    }

    @Nullable
    List<Virus> allViruses();

    @Gson.Ignore
    @Nullable
    default List<Virus> reportableViruses()
    {
        return filterReported(allViruses());
    }

    @Nullable
    PredictedTumorOrigin predictedTumorOrigin();

    private static <T extends Driver> List<T> filterReported(List<T> drivers)
    {
        if (drivers == null) { return null; }
        return drivers.stream().filter(Driver::isReportable).toList();
    }
}
