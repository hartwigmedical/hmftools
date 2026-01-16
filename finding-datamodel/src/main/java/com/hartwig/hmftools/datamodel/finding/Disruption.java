package com.hartwig.hmftools.datamodel.finding;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;

import org.jspecify.annotations.Nullable;

import io.soabase.recordbuilder.core.RecordBuilder;
import jakarta.validation.constraints.NotNull;

@RecordBuilder
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
        @Nullable Breakend breakendStart,
        @Nullable Breakend breakendEnd
) implements Driver
{
    public enum Type
    {
        SOMATIC_DISRUPTION,
        SOMATIC_HOM_DUP_DISRUPTION,
        SOMATIC_HOM_DEL_DISRUPTION,
        GERMLINE_DISRUPTION,
        GERMLINE_HOM_DUP_DISRUPTION;

        public boolean isSomatic()
        {
            return !isGermline();
        }

        public boolean isGermline()
        {
            return this == GERMLINE_DISRUPTION || this == GERMLINE_HOM_DUP_DISRUPTION;
        }

        public boolean isHomozygous()
        {
            return this == SOMATIC_HOM_DUP_DISRUPTION ||
                    this == SOMATIC_HOM_DEL_DISRUPTION ||
                    this == GERMLINE_HOM_DUP_DISRUPTION;
        }
    }

    @NotNull
    @Override
    public String findingKey()
    {
        return driver.findingKey();
    }

    @NotNull
    @Override
    public DriverSource driverSource()
    {
        return driver.driverSource();
    }

    @NotNull
    @Override
    public ReportedStatus reportedStatus()
    {
        return driver.reportedStatus();
    }

    @NotNull
    @Override
    public DriverInterpretation driverInterpretation()
    {
        return driver.driverInterpretation();
    }

    public boolean isHomozygous()
    {
        return type().isHomozygous();
    }

    @NotNull
    public String disruptedRange()
    {
        if(breakendStart != null && breakendEnd != null)
        {
            return exonDescription(breakendStart.exonUp(), breakendStart.exonDown()) + " -> " + exonDescription(breakendEnd.exonUp(),
                    breakendEnd.exonDown());
        }
        if(breakendEnd == null)
        {
            assert breakendStart != null;
            return exonDescription(breakendStart.exonUp(), breakendStart.exonDown()) + " Upstream";
        }
        else
        {
            return exonDescription(breakendEnd.exonUp(), breakendEnd.exonDown()) + " Downstream";
        }
    }

    @NotNull
    private static String exonDescription(int exonUp, int exonDown)
    {
        if(exonUp > 0)
        {
            if(exonUp == exonDown)
            {
                return String.format("Exon %d", exonUp);
            }
            else if(exonDown - exonUp == 1)
            {
                return String.format("Intron %d", exonUp);
            }
        }
        else if(exonUp == 0 && (exonDown == 1 || exonDown == 2))
        {
            return "Promoter Region";
        }
        return String.format("ERROR up=%d, down=%d", exonUp, exonDown);
    }
}
