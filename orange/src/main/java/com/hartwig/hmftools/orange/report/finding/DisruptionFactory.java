package com.hartwig.hmftools.orange.report.finding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.finding.Doubles;
import com.hartwig.hmftools.datamodel.finding.DriverInterpretation;
import com.hartwig.hmftools.datamodel.finding.ImmutableDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class DisruptionFactory {
    private static final String BREAKEND_ORIENTATION_UPSTREAM = "Upstream";
    private static final Logger LOGGER = LogManager.getLogger(DisruptionFactory.class);

    private DisruptionFactory() {
    }

    @NotNull
    public static List<Disruption> convert(@NotNull Collection<LinxBreakend> breakends,
            @NotNull Iterable<LinxSvAnnotation> structuralVariants, boolean hasReliablePurity) {
        List<Disruption> reportableDisruptions = new ArrayList<>();
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = mapBreakendsPerStructuralVariant(breakends);

        for (Pair<LinxBreakend, LinxBreakend> pairedBreakend : pairedMap.values()) {
            LinxBreakend primaryBreakendLeft = pairedBreakend.getLeft();
            LinxBreakend primaryBreakendRight = pairedBreakend.getRight();

            double undisruptedCopyNumber;
            if (primaryBreakendRight != null) {
                undisruptedCopyNumber = Math.min(primaryBreakendLeft.undisruptedCopyNumber(), primaryBreakendRight.undisruptedCopyNumber());

                double copyNumberLeft = primaryBreakendLeft.junctionCopyNumber();
                double copyNumberRight = primaryBreakendRight.junctionCopyNumber();
                if (!Doubles.equal(copyNumberLeft, copyNumberRight)) {
                    LOGGER.warn("The disrupted copy number of a paired sv is not the same on {}", primaryBreakendLeft.gene());
                }
            } else {
                undisruptedCopyNumber = primaryBreakendLeft.undisruptedCopyNumber();
            }
            reportableDisruptions.add(createDisruption(primaryBreakendLeft,
                    rangeField(pairedBreakend),
                    undisruptedCopyNumber,
                    structuralVariants,
                    hasReliablePurity));
        }
        return reportableDisruptions;
    }

    @NotNull
    public static Disruption createDisruption(LinxBreakend primaryBreakendLeft, String disruptedRange, double undisruptedCopyNumber,
            Iterable<LinxSvAnnotation> structuralVariants, boolean hasReliablePurity) {
        return ImmutableDisruption.builder()
                .findingKey(FindingKeys.findingKey(primaryBreakendLeft))
                .isReportable(primaryBreakendLeft.reported())
                .isCandidate(false)
                .driverInterpretation(DriverInterpretation.HIGH)
                .chromosome(primaryBreakendLeft.chromosome())
                .chromosomeBand(primaryBreakendLeft.chromosomeBand())
                .gene(primaryBreakendLeft.gene())
                .type(primaryBreakendLeft.type())
                .disruptedRange(disruptedRange)
                .disruptedCopies(hasReliablePurity ? primaryBreakendLeft.junctionCopyNumber() : null)
                .undisruptedCopies(hasReliablePurity ? undisruptedCopyNumber : null)
                .clusterId(determineClusterId(structuralVariants, primaryBreakendLeft))
                .build();
    }

    @VisibleForTesting
    @Nullable
    static Integer determineClusterId(@NotNull Iterable<LinxSvAnnotation> structuralVariants, @NotNull LinxBreakend breakend) {
        for (LinxSvAnnotation structuralVariant : structuralVariants) {
            if (structuralVariant.svId() == breakend.svId()) {
                return structuralVariant.clusterId();
            }
        }

        LOGGER.warn("Could not find cluster ID for breakend with svId {}", breakend.svId());
        return null;
    }

    @NotNull
    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> mapBreakendsPerStructuralVariant(
            @NotNull Collection<LinxBreakend> breakends) {
        Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript = breakends.stream()
                .collect(Collectors.groupingBy(breakend -> new SvAndTranscriptKey(breakend.svId(), breakend.transcript()),
                        HashMap::new,
                        Collectors.toList()));
        return toPairedMap(breakendsPerSvAndTranscript);
    }

    @NotNull
    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> toPairedMap(
            @NotNull Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript) {
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = new HashMap<>();

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
        return breakend.geneOrientation().equals(BREAKEND_ORIENTATION_UPSTREAM);
    }

    private record SvAndTranscriptKey(int variantId, @NotNull String transcriptId) {
    }
}
