package com.hartwig.hmftools.finding.datamodel;

import java.util.Set;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

@SuppressWarnings("unused")
@RecordBuilder
public record Disruption(
        @NotNull DriverFields driver,
        @NotNull Type type,
        @NotNull String chromosome,
        @NotNull String chromosomeBand,
        @NotNull String gene,
        @NotNull String transcript,
        boolean isCanonical,
        @NotNull BreakendType breakendType,
        double disruptedCopyNumber,
        double undisruptedCopyNumber,
        @Nullable Integer clusterId,
        @Nullable Breakend breakendUp,
        @Nullable Breakend breakendDown
) implements Driver
{
    public enum Type
    {
        DISRUPTION,
        HOM_DUP_DISRUPTION,
        HOM_DEL_DISRUPTION;

        public boolean isHomozygous()
        {
            return this == HOM_DUP_DISRUPTION ||
                    this == HOM_DEL_DISRUPTION;
        }
    }

    public enum BreakendType
    {
        BND,
        DEL,
        DUP,
        INF,
        INS,
        INV,
        SGL
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

    @Override
    public boolean isReported()
    {
        return driver.isReported();
    }

    @Override
    public Set<String> genes()
    {
        return Set.of(gene());
    }

    @Override
    public double driverLikelihood()
    {
        return driver.driverLikelihood();
    }

    public boolean isHomozygous()
    {
        return type().isHomozygous();
    }

    @NotNull
    public String disruptedRange()
    {
        if(breakendUp != null && breakendDown != null)
        {
            return exonDescription(breakendUp.exonUp(), breakendUp.exonDown()) + " -> " + exonDescription(breakendDown.exonUp(),
                    breakendDown.exonDown());
        }
        if(breakendDown == null)
        {
            assert breakendUp != null;
            return exonDescription(breakendUp.exonUp(), breakendUp.exonDown()) + " Upstream";
        }
        else
        {
            return exonDescription(breakendDown.exonUp(), breakendDown.exonDown()) + " Downstream";
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
