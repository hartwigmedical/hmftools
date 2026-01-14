package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public record Disruption(
        @NotNull DriverFields driver,
        @NotNull Type type,
        @NotNull String chromosome,
        @NotNull String chromosomeBand,
        @NotNull String gene,
        @NotNull String transcript,
        boolean isCanonical,
        @NotNull LinxBreakendType breakendType,
        @Nullable Double disruptedCopies,
        @Nullable Double undisruptedCopies,
        @Nullable Integer clusterId,
        @Nullable LinxBreakend breakendStart,
        @Nullable LinxBreakend breakendEnd
) implements Driver
{
    enum Type
    {
        SOMATIC_DISRUPTION,
        SOMATIC_HOM_DUP_DISRUPTION,
        SOMATIC_HOM_DEL_DISRUPTION,
        GERMLINE_DISRUPTION,
        GERMLINE_HOM_DUP_DISRUPTION;

        public boolean isSomatic() { return !isGermline(); }
        public boolean isGermline() { return this == GERMLINE_DISRUPTION || this == GERMLINE_HOM_DUP_DISRUPTION; }
        public boolean isHomozygous() { return this == SOMATIC_HOM_DUP_DISRUPTION ||
                                               this == SOMATIC_HOM_DEL_DISRUPTION ||
                                               this == GERMLINE_HOM_DUP_DISRUPTION; }
    }

    @NotNull @Override public String findingKey() { return driver.findingKey(); }
    @NotNull @Override public DriverSource driverSource() { return driver.driverSource(); }
    @NotNull @Override public ReportedStatus reportedStatus() { return driver.reportedStatus(); }
    @NotNull @Override public DriverInterpretation driverInterpretation() { return driver.driverInterpretation(); }

    boolean isHomozygous() { return type().isHomozygous(); }
}
