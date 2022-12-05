package com.hartwig.hmftools.orange.report.datamodel;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class BreakendEntryFactory {

    private static final Logger LOGGER = LogManager.getLogger(BreakendEntryFactory.class);

    private BreakendEntryFactory() {
    }

    @NotNull
    public static List<BreakendEntry> create(@NotNull List<LinxBreakend> breakends, @NotNull List<LinxSvAnnotation> variants) {
        List<BreakendEntry> entries = Lists.newArrayList();
        for (LinxBreakend breakend : breakends) {
            entries.add(ImmutableBreakendEntry.builder()
                    .location(breakend.chromosome() + breakend.chrBand())
                    .gene(breakend.gene())
                    .canonical(breakend.canonical())
                    .exonUp(breakend.exonUp())
                    .type(breakend.type())
                    .range(range(breakend))
                    .clusterId(determineClusterId(breakend, variants))
                    .junctionCopyNumber(breakend.junctionCopyNumber())
                    .undisruptedCopyNumber(breakend.undisruptedCopyNumber())
                    .build());
        }
        return entries;
    }

    @NotNull
    @VisibleForTesting
    static String range(@NotNull LinxBreakend breakend) {
        String exonRange = null;
        if (breakend.exonUp() > 0) {
            if (breakend.exonUp() == breakend.exonDown()) {
                exonRange = String.format("Exon %d", breakend.exonUp());
            } else if (breakend.exonDown() - breakend.exonUp() == 1) {
                exonRange = String.format("Intron %d", breakend.exonUp());
            }
        } else if (breakend.exonUp() == 0 && (breakend.exonDown() == 1 || breakend.exonDown() == 2)) {
            exonRange = "Promoter Region";
        }

        if (exonRange == null) {
            LOGGER.warn("Could not format range for breakend: {}", breakend);
            return Strings.EMPTY;
        }

        return exonRange + " " + breakend.geneOrientation();
    }

    private static int determineClusterId(@NotNull LinxBreakend breakend, @NotNull List<LinxSvAnnotation> variants) {
        for (LinxSvAnnotation variant : variants) {
            if (variant.svId() == breakend.svId()) {
                return variant.clusterId();
            }
        }

        throw new IllegalStateException("Could not find structural variant that underlies breakend: " + breakend);
    }
}
