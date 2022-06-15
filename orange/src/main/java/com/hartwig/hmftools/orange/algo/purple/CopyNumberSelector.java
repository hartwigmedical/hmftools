package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;
import com.hartwig.hmftools.common.purple.interpretation.ImmutableGainLoss;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public final class CopyNumberSelector {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberSelector.class);

    private CopyNumberSelector() {
    }

    @NotNull
    public static List<GainLoss> selectNearReportableSomaticGains(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, double ploidy,
            @NotNull List<GainLoss> reportableGainsLosses, @NotNull List<DriverGene> driverGenes) {
        List<GainLoss> nearReportableSomaticGains = Lists.newArrayList();
        Set<String> ampDriverGenes = selectAmpDriverGenes(driverGenes);
        for (GeneCopyNumber geneCopyNumber : allGeneCopyNumbers) {
            if (ampDriverGenes.contains(geneCopyNumber.geneName())) {
                double relativeMinCopyNumber = geneCopyNumber.minCopyNumber() / ploidy;
                double relativeMaxCopyNumber = geneCopyNumber.maxCopyNumber() / ploidy;
                if (relativeMinCopyNumber > 2.5 && relativeMaxCopyNumber < 3) {
                    nearReportableSomaticGains.add(toFullGain(geneCopyNumber));
                }
            }
        }

        // Check in case official amp have changed.
        Set<String> reportableGenes = Sets.newHashSet();
        for (GainLoss reportable : reportableGainsLosses) {
            reportableGenes.add(reportable.gene());
        }

        for (GainLoss gain : nearReportableSomaticGains) {
            if (reportableGenes.contains(gain.gene())) {
                LOGGER.warn("Gene {} is selected to be near-reportable but has already been reported!", gain.gene());
            }
        }

        return nearReportableSomaticGains;
    }

    @NotNull
    public static List<GainLoss> selectInterestingUnreportedGainsLosses(@NotNull List<GainLoss> allGainsLosses,
            @NotNull List<GainLoss> reportableGainsLosses) {
        List<GainLoss> unreportedGainLosses = selectUnreportedGainsLosses(allGainsLosses, reportableGainsLosses);

        List<GainLoss> interestingUnreportedGainsLosses = Lists.newArrayList();
        interestingUnreportedGainsLosses.addAll(selectInterestingGains(unreportedGainLosses));
        interestingUnreportedGainsLosses.addAll(selectInterestingLosses(unreportedGainLosses, reportableGainsLosses));
        return interestingUnreportedGainsLosses;
    }

    @NotNull
    private static Set<String> selectAmpDriverGenes(@NotNull List<DriverGene> driverGenes) {
        Set<String> ampGenes = Sets.newHashSet();
        for (DriverGene driverGene : driverGenes) {
            if (driverGene.reportAmplification()) {
                ampGenes.add(driverGene.gene());
            }
        }
        return ampGenes;
    }

    @NotNull
    private static GainLoss toFullGain(@NotNull GeneCopyNumber geneCopyNumber) {
        return ImmutableGainLoss.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.geneName())
                .transcript(geneCopyNumber.transName())
                .isCanonical(geneCopyNumber.isCanonical())
                .interpretation(CopyNumberInterpretation.FULL_GAIN)
                .minCopies(Math.round(Math.max(0, geneCopyNumber.minCopyNumber())))
                .maxCopies(Math.round(Math.max(0, geneCopyNumber.maxCopyNumber())))
                .build();
    }

    @NotNull
    private static List<GainLoss> selectUnreportedGainsLosses(@NotNull List<GainLoss> allGainsLosses,
            @NotNull List<GainLoss> reportableGainsLosses) {
        List<GainLoss> unreportedGainsLosses = Lists.newArrayList();
        for (GainLoss gainLoss : allGainsLosses) {
            if (!reportableGainsLosses.contains(gainLoss)) {
                unreportedGainsLosses.add(gainLoss);
            }
        }
        return unreportedGainsLosses;
    }

    @NotNull
    private static List<GainLoss> selectInterestingGains(@NotNull List<GainLoss> unreportedGainLosses) {
        List<GainLoss> unreportedFullGains = unreportedGainLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN)
                .collect(Collectors.toList());

        Map<CopyNumberKey, GainLoss> maxGainPerLocation = Maps.newHashMap();
        for (GainLoss gain : unreportedFullGains) {
            CopyNumberKey key = new CopyNumberKey(gain.chromosome(), gain.chromosomeBand());
            GainLoss maxGain = maxGainPerLocation.get(key);
            if (maxGain == null || gain.minCopies() > maxGain.minCopies()) {
                maxGainPerLocation.put(key, gain);
            }
        }

        return maxGainPerLocation.values()
                .stream()
                .sorted((o1, o2) -> (int) (o2.minCopies() - o1.minCopies()))
                .collect(Collectors.toList());
    }

    @NotNull
    private static List<GainLoss> selectInterestingLosses(@NotNull List<GainLoss> unreportedGainsLosses,
            @NotNull List<GainLoss> reportableGainsLosses) {
        List<GainLoss> lossesNoAllosomes = Lists.newArrayList();

        List<GainLoss> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        List<GainLoss> reportableLosses = reportableGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        for (GainLoss loss : unreportedLosses) {
            if (!HumanChromosome.fromString(loss.chromosome()).isAllosome() && !locusPresent(reportableLosses,
                    loss.chromosome(),
                    loss.chromosomeBand())) {
                lossesNoAllosomes.add(loss);
            }
        }

        Map<CopyNumberKey, GainLoss> oneLossPerLocation = Maps.newHashMap();
        for (GainLoss loss : lossesNoAllosomes) {
            CopyNumberKey key = new CopyNumberKey(loss.chromosome(), loss.chromosomeBand());
            if (!oneLossPerLocation.containsKey(key)) {
                oneLossPerLocation.put(key, loss);
            }
        }

        return Lists.newArrayList(oneLossPerLocation.values().iterator());
    }

    private static boolean locusPresent(@NotNull List<GainLoss> gainsLosses, @NotNull String chromosome, @NotNull String chromosomeBand) {
        for (GainLoss gainLoss : gainsLosses) {
            if (gainLoss.chromosome().equals(chromosome) && gainLoss.chromosomeBand().equals(chromosomeBand)) {
                return true;
            }
        }

        return false;
    }
}
