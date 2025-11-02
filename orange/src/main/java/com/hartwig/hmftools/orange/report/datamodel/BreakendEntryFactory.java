package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.List;

import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class BreakendEntryFactory
{
    public static final Logger LOGGER = LogManager.getLogger(BreakendEntryFactory.class);

    @NotNull
    public static List<BreakendEntry> create(@NotNull List<LinxBreakend> breakends, @NotNull List<LinxSvAnnotation> variants,
            @NotNull List<LinxDriver> drivers)
    {
        List<BreakendEntry> entries = new ArrayList<>();
        for(LinxBreakend breakend : breakends)
        {
            entries.add(ImmutableBreakendEntry.builder()
                    .linxBreakend(breakend)
                    .location(breakend.chromosome() + breakend.chromosomeBand())
                    .range(range(breakend))
                    .clusterId(determineClusterId(breakend, variants))
                    .undisruptedCopyNumber(correctUndisruptedCopyNumber(breakend, drivers))
                    .build());
        }
        return entries;
    }

    @NotNull
    static String range(@NotNull LinxBreakend breakend)
    {
        String exonRange = null;
        if(breakend.exonUp() > 0)
        {
            if(breakend.exonUp() == breakend.exonDown())
            {
                exonRange = String.format("Exon %d", breakend.exonUp());
            }
            else if(breakend.exonDown() - breakend.exonUp() == 1)
            {
                exonRange = String.format("Intron %d", breakend.exonUp());
            }
        }
        else if(breakend.exonUp() == 0 && (breakend.exonDown() == 1 || breakend.exonDown() == 2))
        {
            exonRange = "Promoter Region";
        }

        if(exonRange == null)
        {
            LOGGER.warn("Could not format range for breakend: {}", breakend);
            return "";
        }

        return exonRange + " " + breakend.geneOrientation();
    }

    private static int determineClusterId(@NotNull LinxBreakend breakend, @NotNull List<LinxSvAnnotation> variants)
    {
        for(LinxSvAnnotation variant : variants)
        {
            if(variant.svId() == breakend.svId())
            {
                return variant.clusterId();
            }
        }

        throw new IllegalStateException("Could not find structural variant that underlies breakend: " + breakend);
    }

    private static double correctUndisruptedCopyNumber(@NotNull LinxBreakend breakend, @NotNull List<LinxDriver> drivers)
    {
        if(breakend.type() == LinxBreakendType.DUP)
        {
            for(LinxDriver driver : drivers)
            {
                if(driver.gene().equals(breakend.gene()) && driver.type() == LinxDriverType.HOM_DUP_DISRUPTION)
                {
                    return Math.max(0.0, breakend.undisruptedCopyNumber() - breakend.junctionCopyNumber());
                }
            }
        }

        return breakend.undisruptedCopyNumber();
    }
}
