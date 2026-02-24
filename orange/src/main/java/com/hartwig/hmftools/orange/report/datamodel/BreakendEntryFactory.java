package com.hartwig.hmftools.orange.report.datamodel;

import static java.lang.Math.max;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.apache.logging.log4j.util.Strings;

public final class BreakendEntryFactory
{
    public static List<BreakendEntry> create(final List<LinxBreakend> breakends)
    {
        List<BreakendEntry> entries = Lists.newArrayList();

        for(LinxBreakend breakend : breakends)
        {
            double undisruptedCopyNumber = breakend.undisruptedCopyNumber();
            double driverLikelihood = 0;

            if(breakend.type() == LinxBreakendType.DUP)
            {
                // TODO: check if this is still required
                undisruptedCopyNumber = max(0.0, breakend.undisruptedCopyNumber() - breakend.junctionCopyNumber());
            }

            entries.add(ImmutableBreakendEntry.builder()
                .location(breakend.chromosome() + breakend.chromosomeBand())
                .gene(breakend.gene())
                .canonical(breakend.isCanonical())
                .exonUp(breakend.exonUp())
                .type(breakend.type())
                .range(range(breakend))
                .junctionCopyNumber(breakend.junctionCopyNumber())
                .undisruptedCopyNumber(undisruptedCopyNumber)
                .driverLikelihood(driverLikelihood)
                .build());
        }

        return entries;
    }

    @VisibleForTesting
    static String range(final LinxBreakend breakend)
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
            return Strings.EMPTY;
        }

        String geneOrientation = switch (breakend.geneOrientation()) {
            case UPSTREAM -> "Upstream";
            case DOWNSTREAM -> "Downstream";
        };

        return exonRange + " " + geneOrientation;
    }
}
