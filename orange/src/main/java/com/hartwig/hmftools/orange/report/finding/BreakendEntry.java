package com.hartwig.hmftools.orange.report.finding;

import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(passAnnotations = { NotNull.class, Nullable.class })
public interface BreakendEntry
{
    @NotNull LinxBreakend linxBreakend();

    @NotNull
    String location();

    @Value.Derived
    @NotNull default String gene()
    {
        return linxBreakend().gene();
    }

    @Value.Derived
    default ReportedStatus reportedStatus()
    {
        return linxBreakend().reportedStatus();
    }

    @Value.Derived
    default boolean canonical()
    {
        return linxBreakend().isCanonical();
    }

    @Value.Derived
    default int exonUp()
    {
        return linxBreakend().exonUp();
    }

    @Value.Derived
    @NotNull default LinxBreakendType type()
    {
        return linxBreakend().type();
    }

    @NotNull
    String range();

    int clusterId();

    @Value.Derived
    default double junctionCopyNumber()
    {
        return linxBreakend().junctionCopyNumber();
    }

    double undisruptedCopyNumber();

}
