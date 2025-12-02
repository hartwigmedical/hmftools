package com.hartwig.hmftools.datamodel.linx;

import java.util.List;

import com.hartwig.hmftools.datamodel.driver.Driver;
import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.finding.Fusion;

import org.immutables.gson.Gson;
import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Gson.TypeAdapters
@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface LinxRecord
{
    @NotNull
    List<LinxDriver> somaticDrivers();

    @NotNull
    List<LinxSvAnnotation> allSomaticStructuralVariants();

    @Nullable
    List<LinxSvAnnotation> allGermlineStructuralVariants();

    @NotNull
    List<Fusion> allSomaticFusions();

    @NotNull
    default List<Fusion> reportableSomaticFusions()
    {
        return allSomaticFusions().stream().filter(Driver::isReportable).toList();
    }

    @NotNull
    List<LinxBreakend> driverSomaticBreakends();

    @Gson.Ignore
    @NotNull
    List<LinxBreakend> otherSomaticBreakends();

    @Gson.Ignore
    @Nullable
    List<LinxBreakend> driverGermlineBreakends();

    @Nullable
    List<LinxBreakend> otherGermlineBreakends();

    @NotNull
    List<Disruption> driverSomaticDisruptions();

    @Nullable
    List<Disruption> driverGermlineDisruptions();

    @NotNull
    List<LinxHomozygousDisruption> somaticHomozygousDisruptions();

    @Nullable
    List<LinxHomozygousDisruption> germlineHomozygousDisruptions();
}
