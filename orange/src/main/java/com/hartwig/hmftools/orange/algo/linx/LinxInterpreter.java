package com.hartwig.hmftools.orange.algo.linx;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.HomozygousDisruption;
import com.hartwig.hmftools.common.linx.ImmutableHomozygousDisruption;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.commons.lang3.tuple.Pair;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LinxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LinxInterpreter.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public LinxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache) {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
    }

    @NotNull
    public LinxInterpretedData interpret(@NotNull LinxData linx) {
        List<LinxFusion> additionalSuspectSomaticFusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(linx.allSomaticFusions(), driverGenes);
        LOGGER.info(" Found an additional {} suspect fusions that are potentially interesting", additionalSuspectSomaticFusions.size());

        List<LinxBreakend> additionalSuspectSomaticBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allSomaticBreakends(),
                        linx.reportableSomaticFusions(),
                        knownFusionCache);
        LOGGER.info(" Found an additional {} suspect breakends that are potentially interesting", additionalSuspectSomaticBreakends.size());

        List<LinxBreakend> reportableGermlineBreakends = filterForPresenceInTumor(linx.reportableGermlineBreakends());
        if (reportableGermlineBreakends != null) {
            LOGGER.info(" Resolved {} reportable germline breakends from {} disruptive breakends",
                    reportableGermlineBreakends.size(),
                    linx.reportableGermlineBreakends().size());
        }

        List<HomozygousDisruption> germlineHomozygousDisruptions = resolveHomozygousDisruptions(reportableGermlineBreakends);
        if (germlineHomozygousDisruptions != null) {
            LOGGER.info(" Resolved {} germline homozygous disruptions from {} breakends",
                    germlineHomozygousDisruptions.size(),
                    reportableGermlineBreakends.size());
        }

        return ImmutableLinxInterpretedData.builder()
                .allSomaticStructuralVariants(linx.allSomaticStructuralVariants())
                .allSomaticFusions(linx.allSomaticFusions())
                .reportableSomaticFusions(linx.reportableSomaticFusions())
                .additionalSuspectSomaticFusions(additionalSuspectSomaticFusions)
                .allSomaticBreakends(linx.allSomaticBreakends())
                .reportableSomaticBreakends(linx.reportableSomaticBreakends())
                .additionalSuspectSomaticBreakends(additionalSuspectSomaticBreakends)
                .somaticHomozygousDisruptions(linx.somaticHomozygousDisruptions())
                .allGermlineStructuralVariants(linx.allGermlineStructuralVariants())
                .allGermlineBreakends(linx.allGermlineBreakends())
                .reportableGermlineBreakends(reportableGermlineBreakends)
                .germlineHomozygousDisruptions(germlineHomozygousDisruptions)
                .build();
    }

    @Nullable
    @VisibleForTesting
    static List<LinxBreakend> filterForPresenceInTumor(@Nullable List<LinxBreakend> breakends) {
        if (breakends == null) {
            return null;
        }

        List<LinxBreakend> filtered = Lists.newArrayList();
        for (LinxBreakend breakend : breakends) {
            if (breakend.junctionCopyNumber() > 0.1) {
                filtered.add(breakend);
            }
        }

        return filtered;
    }

    @Nullable
    @VisibleForTesting
    static List<HomozygousDisruption> resolveHomozygousDisruptions(@Nullable List<LinxBreakend> breakends) {
        if (breakends == null) {
            return null;
        }

        List<HomozygousDisruption> homozygousDisruptions = Lists.newArrayList();
        for (Pair<LinxBreakend, LinxBreakend> breakendPair : createBreakendPairsPerSvId(breakends)) {
            LinxBreakend first = breakendPair.getLeft();
            LinxBreakend second = breakendPair.getRight();

            boolean bothReported = first.reportedDisruption() && second.reportedDisruption();
            boolean bothDup = first.type() == StructuralVariantType.DUP && second.type() == StructuralVariantType.DUP;
            boolean bothSameGene = first.gene().equals(second.gene());
            boolean bothSameTranscript = first.transcriptId().equals(second.transcriptId());
            boolean bothDisruptRemainingCopies = first.junctionCopyNumber() >= first.undisruptedCopyNumber()
                    && second.junctionCopyNumber() >= second.undisruptedCopyNumber();

            if (bothReported && bothDup && bothSameGene && bothSameTranscript && bothDisruptRemainingCopies) {
                homozygousDisruptions.add(ImmutableHomozygousDisruption.builder()
                        .chromosome(first.chromosome())
                        .chromosomeBand(first.chrBand())
                        .gene(first.gene())
                        .transcript(first.transcriptId())
                        .isCanonical(first.canonical())
                        .build());
            }
        }
        return homozygousDisruptions;
    }

    @NotNull
    private static Collection<Pair<LinxBreakend, LinxBreakend>> createBreakendPairsPerSvId(@NotNull List<LinxBreakend> breakends) {
        Map<Integer, Pair<LinxBreakend, LinxBreakend>> mapPerSvId = Maps.newHashMap();
        for (LinxBreakend breakend : breakends) {
            if (!mapPerSvId.containsKey(breakend.svId())) {
                LinxBreakend paired = findPaired(breakends, breakend);
                if (paired != null) {
                    mapPerSvId.put(breakend.svId(), Pair.of(breakend, paired));
                }
            }
        }
        return mapPerSvId.values();
    }

    @Nullable
    private static LinxBreakend findPaired(@NotNull List<LinxBreakend> breakends, @NotNull LinxBreakend breakendToFindPairFor) {
        for (LinxBreakend breakend : breakends) {
            if (breakend != breakendToFindPairFor && breakend.svId() == breakendToFindPairFor.svId()) {
                return breakend;
            }
        }
        return null;
    }
}
