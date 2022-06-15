package com.hartwig.hmftools.common.linx;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.sv.linx.LinxBreakend;
import com.hartwig.hmftools.common.sv.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ReportableGeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(ReportableGeneDisruptionFactory.class);

    private ReportableGeneDisruptionFactory() {
    }

    @NotNull
    public static List<ReportableGeneDisruption> convert(@NotNull List<LinxBreakend> breakends,
            @NotNull List<LinxSvAnnotation> structuralVariants) {
        List<ReportableGeneDisruption> reportableDisruptions = Lists.newArrayList();
        Map<SvAndGeneKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = mapBreakendsPerStructuralVariant(breakends);

        for (Pair<LinxBreakend, LinxBreakend> pairedBreakend : pairedMap.values()) {
            LinxBreakend primaryBreakendLeft = pairedBreakend.getLeft();
            LinxBreakend primaryBreakendRight = pairedBreakend.getRight();

            if (primaryBreakendRight != null) {
                double lowestUndisruptedCopyNumber =
                        Math.min(primaryBreakendLeft.undisruptedCopyNumber(), primaryBreakendRight.undisruptedCopyNumber());

                Double copyNumberLeft = primaryBreakendLeft.junctionCopyNumber();
                Double copyNumberRight = primaryBreakendRight.junctionCopyNumber();
                if (copyNumberLeft != null && copyNumberRight != null && !Doubles.equal(copyNumberLeft, copyNumberRight)) {
                    LOGGER.warn("The disrupted copy number of a paired sv is not the same on {}", primaryBreakendLeft.gene());
                }
                reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                        .location(primaryBreakendLeft.chromosome() + primaryBreakendLeft.chrBand())
                        .gene(primaryBreakendLeft.gene())
                        .transcriptId(primaryBreakendLeft.transcriptId())
                        .isCanonical(primaryBreakendLeft.canonical())
                        .type(primaryBreakendLeft.type())
                        .range(rangeField(pairedBreakend))
                        .junctionCopyNumber(primaryBreakendLeft.junctionCopyNumber())
                        .undisruptedCopyNumber(Math.max(0, lowestUndisruptedCopyNumber))
                        .firstAffectedExon(primaryBreakendLeft.exonUp())
                        .clusterId(determineClusterId(structuralVariants, primaryBreakendLeft))
                        .build());
            } else {
                reportableDisruptions.add(ImmutableReportableGeneDisruption.builder()
                        .location(primaryBreakendLeft.chromosome() + primaryBreakendLeft.chrBand())
                        .gene(primaryBreakendLeft.gene())
                        .transcriptId(primaryBreakendLeft.transcriptId())
                        .isCanonical(primaryBreakendLeft.canonical())
                        .type(primaryBreakendLeft.type())
                        .range(rangeField(pairedBreakend))
                        .junctionCopyNumber(primaryBreakendLeft.junctionCopyNumber())
                        .undisruptedCopyNumber(Math.max(0, primaryBreakendLeft.undisruptedCopyNumber()))
                        .firstAffectedExon(primaryBreakendLeft.exonUp())
                        .clusterId(determineClusterId(structuralVariants, primaryBreakendLeft))
                        .build());
            }
        }
        return reportableDisruptions;
    }

    @VisibleForTesting
    @Nullable
    static Integer determineClusterId(@NotNull List<LinxSvAnnotation> structuralVariants, @NotNull LinxBreakend breakend) {
        for (LinxSvAnnotation structuralVariant : structuralVariants) {
            if (structuralVariant.svId() == breakend.svId()) {
                return structuralVariant.clusterId();
            }
        }

        if (!structuralVariants.isEmpty()) {
            LOGGER.warn("Could not find cluster ID for breakend with svId {}", breakend.svId());
        }

        return null;
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<LinxBreakend, LinxBreakend>> mapBreakendsPerStructuralVariant(
            @NotNull List<LinxBreakend> breakends) {
        Map<SvAndGeneKey, List<LinxBreakend>> breakendsPerSvAndGene = Maps.newHashMap();
        for (LinxBreakend breakend : breakends) {
            SvAndGeneKey key = new SvAndGeneKey(breakend.svId(), breakend.gene());
            List<LinxBreakend> currentBreakends = breakendsPerSvAndGene.get(key);
            if (currentBreakends == null) {
                currentBreakends = Lists.newArrayList();
            }
            currentBreakends.add(breakend);
            breakendsPerSvAndGene.put(key, currentBreakends);
        }

        return toPairedMap(breakendsPerSvAndGene);
    }

    @NotNull
    private static Map<SvAndGeneKey, Pair<LinxBreakend, LinxBreakend>> toPairedMap(
            @NotNull Map<SvAndGeneKey, List<LinxBreakend>> breakendsPerSvAndGene) {
        Map<SvAndGeneKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = Maps.newHashMap();

        for (Map.Entry<SvAndGeneKey, List<LinxBreakend>> entry : breakendsPerSvAndGene.entrySet()) {
            List<LinxBreakend> breakends = entry.getValue();

            if (breakends.size() != 1 && breakends.size() != 2) {
                LOGGER.warn("Found unusual number of breakends on single event: {}", breakends.size());
                continue;
            }

            LinxBreakend left;
            LinxBreakend right;
            if (breakends.size() == 1) {
                left = breakends.get(0);
                right = null;
            } else {
                boolean firstBeforeSecond = breakends.get(0).exonUp() <= breakends.get(1).exonUp();
                left = firstBeforeSecond ? breakends.get(0) : breakends.get(1);
                right = firstBeforeSecond ? breakends.get(1) : breakends.get(0);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static String rangeField(@NotNull Pair<LinxBreakend, LinxBreakend> pairedBreakend) {
        LinxBreakend primary = pairedBreakend.getLeft();
        LinxBreakend secondary = pairedBreakend.getRight();

        if (secondary == null) {
            return exonDescription(primary.exonUp(), primary.exonDown()) + (isUpstream(primary) ? " Upstream" : " Downstream");
        } else {
            return exonDescription(primary.exonUp(), primary.exonDown()) + " -> " + exonDescription(secondary.exonUp(),
                    secondary.exonDown());
        }
    }

    @NotNull
    private static String exonDescription(int exonUp, int exonDown) {
        if (exonUp > 0) {
            if (exonUp == exonDown) {
                return String.format("Exon %d", exonUp);
            } else if (exonDown - exonUp == 1) {
                return String.format("Intron %d", exonUp);
            }
        } else if (exonUp == 0 && (exonDown == 1 || exonDown == 2)) {
            return "Promoter Region";
        }

        return String.format("ERROR up=%d, down=%d", exonUp, exonDown);
    }

    private static boolean isUpstream(@NotNull LinxBreakend breakend) {
        return breakend.orientation() * breakend.strand() < 0;
    }

    private static class SvAndGeneKey {
        private final int variantId;
        @NotNull
        private final String gene;

        private SvAndGeneKey(final int variantId, @NotNull final String gene) {
            this.variantId = variantId;
            this.gene = gene;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SvAndGeneKey that = (SvAndGeneKey) o;
            return Objects.equals(variantId, that.variantId) && Objects.equals(gene, that.gene);
        }

        @Override
        public int hashCode() {
            return Objects.hash(variantId, gene);
        }
    }
}
