package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.Map;
import java.util.Objects;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.utils.Doubles;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class GeneDisruptionFactory {

    private static final Logger LOGGER = LogManager.getLogger(GeneDisruptionFactory.class);

    private GeneDisruptionFactory() {
    }

    @NotNull
    public static List<GeneDisruption> convert(final List<LinxBreakend> breakends, final List<LinxSvAnnotation> structuralVariants) {
        List<GeneDisruption> reportableDisruptions = Lists.newArrayList();
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = mapBreakendsPerStructuralVariant(breakends);

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
                reportableDisruptions.add(ImmutableGeneDisruption.builder()
                        .location(primaryBreakendLeft.chromosome() + primaryBreakendLeft.chrBand())
                        .gene(primaryBreakendLeft.gene())
                        .transcriptId(primaryBreakendLeft.transcriptId())
                        .isCanonical(primaryBreakendLeft.canonical())
                        .type(primaryBreakendLeft.type().toString())
                        .range(rangeField(pairedBreakend))
                        .junctionCopyNumber(primaryBreakendLeft.junctionCopyNumber())
                        .undisruptedCopyNumber(Math.max(0, lowestUndisruptedCopyNumber))
                        .firstAffectedExon(primaryBreakendLeft.exonUp())
                        .clusterId(determineClusterId(structuralVariants, primaryBreakendLeft))
                        .build());
            } else {
                reportableDisruptions.add(ImmutableGeneDisruption.builder()
                        .location(primaryBreakendLeft.chromosome() + primaryBreakendLeft.chrBand())
                        .gene(primaryBreakendLeft.gene())
                        .transcriptId(primaryBreakendLeft.transcriptId())
                        .isCanonical(primaryBreakendLeft.canonical())
                        .type(primaryBreakendLeft.type().toString())
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
    static Integer determineClusterId(final List<LinxSvAnnotation> structuralVariants, final LinxBreakend breakend) {
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

    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> mapBreakendsPerStructuralVariant(
            final List<LinxBreakend> breakends) {
        Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript = Maps.newHashMap();
        for (LinxBreakend breakend : breakends) {
            SvAndTranscriptKey key = new SvAndTranscriptKey(breakend.svId(), breakend.transcriptId());
            List<LinxBreakend> currentBreakends = breakendsPerSvAndTranscript.get(key);
            if (currentBreakends == null) {
                currentBreakends = Lists.newArrayList();
            }
            currentBreakends.add(breakend);
            breakendsPerSvAndTranscript.put(key, currentBreakends);
        }

        return toPairedMap(breakendsPerSvAndTranscript);
    }

    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> toPairedMap(
            final Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript) {
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = Maps.newHashMap();

        for (Map.Entry<SvAndTranscriptKey, List<LinxBreakend>> entry : breakendsPerSvAndTranscript.entrySet()) {
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

    private static String rangeField(final Pair<LinxBreakend, LinxBreakend> pairedBreakend) {
        LinxBreakend primary = pairedBreakend.getLeft();
        LinxBreakend secondary = pairedBreakend.getRight();

        if (secondary == null) {
            return exonDescription(primary.exonUp(), primary.exonDown()) + (isUpstream(primary) ? " Upstream" : " Downstream");
        } else {
            return exonDescription(primary.exonUp(), primary.exonDown()) + " -> " + exonDescription(secondary.exonUp(),
                    secondary.exonDown());
        }
    }

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

    private static boolean isUpstream(final LinxBreakend breakend) {
        return breakend.orientation() * breakend.strand() < 0;
    }

    private static class SvAndTranscriptKey {
        private final int variantId;
        @NotNull
        private final String transcriptId;

        private SvAndTranscriptKey(final int variantId, @NotNull final String transcriptId) {
            this.variantId = variantId;
            this.transcriptId = transcriptId;
        }

        @Override
        public boolean equals(final Object o) {
            if (this == o) {
                return true;
            }
            if (o == null || getClass() != o.getClass()) {
                return false;
            }
            final SvAndTranscriptKey that = (SvAndTranscriptKey) o;
            return Objects.equals(variantId, that.variantId) && Objects.equals(transcriptId, that.transcriptId);
        }

        @Override
        public int hashCode() {
            return Objects.hash(variantId, transcriptId);
        }
    }
}
