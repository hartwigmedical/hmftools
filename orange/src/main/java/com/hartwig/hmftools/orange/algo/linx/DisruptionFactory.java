package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.finding.Disruption;
import com.hartwig.hmftools.datamodel.finding.Doubles;
import com.hartwig.hmftools.datamodel.finding.ImmutableDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;
import com.hartwig.hmftools.orange.algo.util.FindingKeys;

import org.apache.commons.lang3.tuple.Pair;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// create Disruption findings from Linx breakends
public class DisruptionFactory
{
    @NotNull
    public static List<Disruption> createDisruptions(
            @NotNull FindingKeys.SampleType sampleType,
            @NotNull Collection<LinxBreakend> breakends,
            @NotNull Collection<LinxSvAnnotation> structuralVariants,
            boolean hasReliablePurity)
    {
        List<Disruption> reportableDisruptions = new ArrayList<>();
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = mapBreakendsPerStructuralVariant(breakends);

        for(Pair<LinxBreakend, LinxBreakend> pairedBreakend : pairedMap.values())
        {
            LinxBreakend primarybreakendStart = pairedBreakend.getLeft();
            LinxBreakend primarybreakendEnd = pairedBreakend.getRight();

            double undisruptedCopyNumber;
            if(primarybreakendStart != null && primarybreakendEnd != null)
            {
                undisruptedCopyNumber = Math.min(primarybreakendStart.undisruptedCopyNumber(), primarybreakendEnd.undisruptedCopyNumber());

                double copyNumberLeft = primarybreakendStart.junctionCopyNumber();
                double copyNumberRight = primarybreakendEnd.junctionCopyNumber();
                if (!Doubles.equal(copyNumberLeft, copyNumberRight))
                {
                    LOGGER.warn("The disrupted copy number of a paired sv is not the same on {}", primarybreakendStart.gene());
                }
            }
            else
            {
                LinxBreakend breakend = primarybreakendStart == null ? primarybreakendEnd : primarybreakendStart;
                undisruptedCopyNumber = breakend.undisruptedCopyNumber();
            }
            reportableDisruptions.add(createDisruption(
                    sampleType,
                    primarybreakendStart,
                    primarybreakendEnd,
                    undisruptedCopyNumber,
                    structuralVariants,
                    hasReliablePurity));
        }
        return reportableDisruptions;
    }

    @NotNull
    public static Disruption createDisruption(
            @NotNull FindingKeys.SampleType sampleType,
            @Nullable LinxBreakend breakendStart,
            @Nullable LinxBreakend breakendEnd,
            double undisruptedCopyNumber,
            Collection<LinxSvAnnotation> structuralVariants, boolean hasReliablePurity)
    {
        LinxBreakend breakend = breakendStart == null ? breakendEnd : breakendStart;

        if(breakend == null)
        {
            // should not be possible
            throw new IllegalStateException("Disruption with no breakend");
        }

        ReportedStatus reportedStatus = breakend.reportedStatus();

        return ImmutableDisruption.builder()
                .findingKey(FindingKeys.disruption(sampleType, breakend))
                .reportedStatus(reportedStatus)
                .driverInterpretation(DriverInterpretation.HIGH) // TODOHWL: fix
                .chromosome(breakend.chromosome())
                .chromosomeBand(breakend.chromosomeBand())
                .gene(breakend.gene())
                .transcript(breakend.transcript())
                .type(breakend.type())
                .disruptedCopies(hasReliablePurity ? breakend.junctionCopyNumber() : null)
                .undisruptedCopies(hasReliablePurity ? undisruptedCopyNumber : null)
                .clusterId(determineClusterId(structuralVariants, breakend))
                .breakendStart(breakendStart)
                .breakendEnd(breakendEnd)
                .build();
    }

    @VisibleForTesting
    @Nullable
    static Integer determineClusterId(@NotNull Collection<LinxSvAnnotation> structuralVariants, @NotNull LinxBreakend breakend)
    {
        Optional<LinxSvAnnotation> sv = structuralVariants.stream().filter(o -> o.svId() == breakend.svId()).findFirst();
        if (sv.isPresent())
        {
            return sv.get().clusterId();
        }
        else
        {
            LOGGER.warn("Could not find cluster ID for breakend with svId {}", breakend.svId());
            return null;
        }
    }

    @NotNull
    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> mapBreakendsPerStructuralVariant(
            @NotNull Collection<LinxBreakend> breakends)
    {
        Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript = breakends.stream()
                .collect(Collectors.groupingBy(breakend -> new SvAndTranscriptKey(breakend.svId(), breakend.transcript()),
                        HashMap::new,
                        Collectors.toList()));
        return toPairedMap(breakendsPerSvAndTranscript);
    }

    @NotNull
    private static Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> toPairedMap(
            @NotNull Map<SvAndTranscriptKey, List<LinxBreakend>> breakendsPerSvAndTranscript)
    {
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = new HashMap<>();

        for (Map.Entry<SvAndTranscriptKey, List<LinxBreakend>> entry : breakendsPerSvAndTranscript.entrySet())
        {
            List<LinxBreakend> breakends = entry.getValue();

            if (breakends.size() != 1 && breakends.size() != 2)
            {
                LOGGER.warn("Found unusual number of breakends on single event: {}", breakends.size());
                continue;
            }

            LinxBreakend left;
            LinxBreakend right;
            if (breakends.size() == 1)
            {
                LinxBreakend breakend = breakends.get(0);
                if(breakend.geneOrientation() == LinxGeneOrientation.UPSTREAM)
                {
                    left = breakend;
                    right = null;
                }
                else
                {
                    right = breakend;
                    left = null;
                }
            }
            else
            {
                boolean firstBeforeSecond = breakends.get(0).exonUp() <= breakends.get(1).exonUp();
                left = firstBeforeSecond ? breakends.get(0) : breakends.get(1);
                right = firstBeforeSecond ? breakends.get(1) : breakends.get(0);
            }

            pairedMap.put(entry.getKey(), Pair.of(left, right));
        }

        return pairedMap;
    }

    @NotNull
    private static String rangeField(@Nullable LinxBreakend breakendStart, @Nullable LinxBreakend breakendEnd)
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
        if (exonUp > 0)
        {
            if (exonUp == exonDown)
            {
                return String.format("Exon %d", exonUp);
            }
            else if (exonDown - exonUp == 1)
            {
                return String.format("Intron %d", exonUp);
            }
        }
        else if (exonUp == 0 && (exonDown == 1 || exonDown == 2))
        {
            return "Promoter Region";
        }

        return String.format("ERROR up=%d, down=%d", exonUp, exonDown);
    }

    private static boolean isUpstream(@NotNull LinxBreakend breakend)
    {
        return breakend.geneOrientation() == LinxGeneOrientation.UPSTREAM;
    }

    private record SvAndTranscriptKey(int variantId, @NotNull String transcriptId)
    {
    }
}
