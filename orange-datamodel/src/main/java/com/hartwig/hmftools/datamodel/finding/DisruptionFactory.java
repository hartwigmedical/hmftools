package com.hartwig.hmftools.datamodel.finding;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;

import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.DriverSource;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxDriver;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

// create Disruption findings from Linx breakends
public class DisruptionFactory
{
    // when can we move to kotlin
    @FunctionalInterface
    private interface DisruptionTypeFinder {
        Disruption.Type apply(String gene, String transcript);
    }

    public static final Logger LOGGER = LogManager.getLogger(DisruptionFactory.class);

    private record Pair<A, B>(@Nullable A left, @Nullable B right) {}

    // for the types of disruptions, see hmftools.common.driver.DriverType, they are:
    // germline: GERMLINE_DISRUPTION, GERMLINE_HOM_DUP_DISRUPTION
    // somatic: HOM_DUP_DISRUPTION, HOM_DEL_DISRUPTION, DISRUPTION

    // we can work out which disruption type it was by the following logic:
    // Somatic:
    //  1. check if the disruption exists in the drivers list, and if it is, use the driver type if the type is HOM_DUP_DISRUPTION or HOM_DEL_DISRUPTION.
    //  2. if not, use SOMATIC_DISRUPTION as the type.
    // Germline:
    //  1. check if the disruption exists in the homozygous disruption list, if it is, then type is GERMLINE_HOM_DUP_DISRUPTION.
    //  2. if not, then type is GERMLINE_DISRUPTION.

    // more note on Homozygous disruptions:
    // right now it is created from linx driver catalog entries instead of the linx breakends.
    // It select HOM_DUP_DISRUPTION, HOM_DEL_DISRUPTION and GERMLINE_HOM_DUP_DISRUPTION. Unfortunately there is no foolproof way to link
    // back to the breakends, it probably will work by just selecting first reportable disruption with the same gene / transcript.
    // We can do it for the backport version if that makes it easier.

    @NotNull
    public static List<Disruption> createGermlineDisruptions(
            @NotNull Collection<LinxBreakend> breakends,
            @NotNull Collection<LinxSvAnnotation> structuralVariants,
            @NotNull List<LinxHomozygousDisruption> germlineHomozygousDisruptions)
    {
        DisruptionTypeFinder findDisruptionType = (gene, transcript) ->
            germlineHomozygousDisruptions.stream()
                .anyMatch(o -> o.gene().equals(gene) && o.transcript().equals(transcript))
                ? Disruption.Type.GERMLINE_HOM_DUP_DISRUPTION
                : Disruption.Type.GERMLINE_DISRUPTION;

        return createDisruptions(DriverSource.GERMLINE, breakends, structuralVariants, true, findDisruptionType);
    }

    @NotNull
    public static List<Disruption> createSomaticDisruptions(
            @NotNull Collection<LinxBreakend> breakends,
            @NotNull Collection<LinxSvAnnotation> structuralVariants,
            @NotNull List<LinxDriver> linxDrivers,
            boolean hasReliablePurity)
    {
        Map<String, Disruption.Type> geneDriverTypeMap = new HashMap<>();
        for(LinxDriver linxDriver : linxDrivers)
        {
            Disruption.Type disruptionType = switch (linxDriver.type())
            {
                case HOM_DUP_DISRUPTION -> Disruption.Type.SOMATIC_HOM_DUP_DISRUPTION;
                case HOM_DEL_DISRUPTION -> Disruption.Type.SOMATIC_HOM_DEL_DISRUPTION;
                default -> Disruption.Type.SOMATIC_DISRUPTION;
            };
            geneDriverTypeMap.put(linxDriver.gene(), disruptionType);
        }
        DisruptionTypeFinder findDisruptionType = (gene, transcript) ->
            geneDriverTypeMap.getOrDefault(gene, Disruption.Type.SOMATIC_DISRUPTION);

        return createDisruptions(DriverSource.SOMATIC, breakends, structuralVariants, hasReliablePurity, findDisruptionType);
    }

    @NotNull
    private static List<Disruption> createDisruptions(
                @NotNull DriverSource sampleType,
                @NotNull Collection<LinxBreakend> breakends,
                @NotNull Collection<LinxSvAnnotation> structuralVariants,
                boolean hasReliablePurity,
                DisruptionTypeFinder disruptionTypeFinder)
    {
        List<Disruption> reportableDisruptions = new ArrayList<>();
        Map<SvAndTranscriptKey, Pair<LinxBreakend, LinxBreakend>> pairedMap = mapBreakendsPerStructuralVariant(breakends);

        for(Pair<LinxBreakend, LinxBreakend> pairedBreakend : pairedMap.values())
        {
            LinxBreakend primaryBreakendStart = pairedBreakend.left();
            LinxBreakend primaryBreakendEnd = pairedBreakend.right();

            LinxBreakend breakend = primaryBreakendStart == null ? primaryBreakendEnd : primaryBreakendStart;
            assert breakend != null;

            double undisruptedCopyNumber;
            if(primaryBreakendStart != null && primaryBreakendEnd != null)
            {
                undisruptedCopyNumber = Math.min(primaryBreakendStart.undisruptedCopyNumber(), primaryBreakendEnd.undisruptedCopyNumber());

                double copyNumberLeft = primaryBreakendStart.junctionCopyNumber();
                double copyNumberRight = primaryBreakendEnd.junctionCopyNumber();
                if (!Doubles.equal(copyNumberLeft, copyNumberRight))
                {
                    LOGGER.warn("The disrupted copy number of a paired sv is not the same on {}", primaryBreakendStart.gene());
                }
            }
            else
            {
                undisruptedCopyNumber = breakend.undisruptedCopyNumber();
            }

            Disruption.Type disruptionType = disruptionTypeFinder.apply(breakend.gene(), breakend.transcript());

            reportableDisruptions.add(createDisruption(
                    sampleType,
                    disruptionType,
                    primaryBreakendStart,
                    primaryBreakendEnd,
                    undisruptedCopyNumber,
                    structuralVariants,
                    hasReliablePurity));
        }
        return reportableDisruptions;
    }

    @NotNull
    public static Disruption createDisruption(
            @NotNull DriverSource sourceSample,
            @NotNull Disruption.Type disruptionType,
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

        ReportedStatus reportedStatus = breakend.reported() ? ReportedStatus.REPORTED : ReportedStatus.NOT_REPORTED;

        return ImmutableDisruption.builder()
                .findingKey(FindingKeys.disruption(sourceSample, breakend))
                .driverSource(sourceSample)
                .reportedStatus(reportedStatus)
                .driverInterpretation(DriverInterpretation.HIGH) // TODOHWL: fix
                .type(disruptionType)
                .chromosome(breakend.chromosome())
                .chromosomeBand(breakend.chromosomeBand())
                .gene(breakend.gene())
                .isCanonical(breakend.isCanonical())
                .transcript(breakend.transcript())
                .breakendType(breakend.type())
                .disruptedCopies(hasReliablePurity ? breakend.junctionCopyNumber() : null)
                .undisruptedCopies(hasReliablePurity ? undisruptedCopyNumber : null)
                .clusterId(determineClusterId(structuralVariants, breakend))
                .breakendStart(breakendStart)
                .breakendEnd(breakendEnd)
                .build();
    }

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
                if(breakend.geneOrientation() == LinxGeneOrientation.Upstream)
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

            pairedMap.put(entry.getKey(), new Pair<>(left, right));
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
        return breakend.geneOrientation() == LinxGeneOrientation.Upstream;
    }

    private record SvAndTranscriptKey(int variantId, @NotNull String transcriptId)
    {
    }
}
