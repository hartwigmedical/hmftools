package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

final class CopyNumberSelector {

    private static final Logger LOGGER = LogManager.getLogger(CopyNumberSelector.class);

    private CopyNumberSelector() {
    }

    @NotNull
    public static List<PurpleGainLoss> selectNearReportableSomaticGains(@NotNull List<GeneCopyNumber> allGeneCopyNumbers, double ploidy,
            @NotNull List<PurpleGainLoss> reportableGainsLosses, @NotNull List<DriverGene> driverGenes) {
        List<PurpleGainLoss> nearReportableSomaticGains = Lists.newArrayList();
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
        for (PurpleGainLoss reportable : reportableGainsLosses) {
            reportableGenes.add(reportable.gene());
        }

        for (PurpleGainLoss gain : nearReportableSomaticGains) {
            if (reportableGenes.contains(gain.gene())) {
                LOGGER.warn("Gene {} is selected to be near-reportable but has already been reported!", gain.gene());
            }
        }

        return nearReportableSomaticGains;
    }

    @NotNull
    public static List<PurpleGainLoss> selectInterestingUnreportedGainsLosses(@NotNull List<PurpleGainLoss> allGainsLosses,
            @NotNull List<PurpleGainLoss> reportableGainsLosses) {
        List<PurpleGainLoss> unreportedGainLosses = selectUnreportedGainsLosses(allGainsLosses, reportableGainsLosses);

        List<PurpleGainLoss> interestingUnreportedGainsLosses = Lists.newArrayList();
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
    private static PurpleGainLoss toFullGain(@NotNull GeneCopyNumber geneCopyNumber) {
        return ImmutablePurpleGainLoss.builder()
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
    private static List<PurpleGainLoss> selectUnreportedGainsLosses(@NotNull List<PurpleGainLoss> allGainsLosses,
            @NotNull List<PurpleGainLoss> reportableGainsLosses) {
        List<PurpleGainLoss> unreportedGainsLosses = Lists.newArrayList();
        for (PurpleGainLoss gainLoss : allGainsLosses) {
            if (!reportableGainsLosses.contains(gainLoss)) {
                unreportedGainsLosses.add(gainLoss);
            }
        }
        return unreportedGainsLosses;
    }

    @NotNull
    private static List<PurpleGainLoss> selectInterestingGains(@NotNull List<PurpleGainLoss> unreportedGainLosses) {
        List<PurpleGainLoss> unreportedFullGains = unreportedGainLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN)
                .collect(Collectors.toList());

        Map<CopyNumberKey, PurpleGainLoss> bestGainPerLocation = Maps.newHashMap();
        for (PurpleGainLoss gain : unreportedFullGains) {
            CopyNumberKey key = new CopyNumberKey(gain.chromosome(), gain.chromosomeBand());
            PurpleGainLoss bestGain = bestGainPerLocation.get(key);
            if (bestGain == null) {
                bestGainPerLocation.put(key, gain);
            } else {
                if (gain.minCopies() > bestGain.minCopies()) {
                    bestGainPerLocation.put(key, gain);
                }
            }
        }

        return Lists.newArrayList(bestGainPerLocation.values().iterator());
    }

    @NotNull
    private static List<PurpleGainLoss> selectInterestingLosses(@NotNull List<PurpleGainLoss> unreportedGainsLosses,
            @NotNull List<PurpleGainLoss> reportableGainsLosses) {
        List<PurpleGainLoss> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        List<PurpleGainLoss> reportableLosses = reportableGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        List<PurpleGainLoss> lossesAutosomes = Lists.newArrayList();
        for (PurpleGainLoss loss : unreportedLosses) {
            if (HumanChromosome.fromString(loss.chromosome()).isAutosome()) {
                if (!locusPresent(reportableLosses, loss.chromosome(), loss.chromosomeBand())) {
                    lossesAutosomes.add(loss);
                }
            }
        }

        Map<CopyNumberKey, PurpleGainLoss> bestLossPerLocation = Maps.newHashMap();
        for (PurpleGainLoss loss : lossesAutosomes) {
            CopyNumberKey key = new CopyNumberKey(loss.chromosome(), loss.chromosomeBand());
            PurpleGainLoss bestLoss = bestLossPerLocation.get(key);
            if (bestLoss == null) {
                bestLossPerLocation.put(key, loss);
            } else {
                boolean pickOtherWhenEqual = bestLoss.gene().compareTo(loss.gene()) <= 0;
                if (pickOtherWhenEqual) {
                    bestLossPerLocation.put(key, loss);
                }
            }
        }

        return Lists.newArrayList(bestLossPerLocation.values().iterator());
    }

    private static boolean locusPresent(@NotNull List<PurpleGainLoss> gainsLosses, @NotNull String chromosome, @NotNull String chromosomeBand) {
        for (PurpleGainLoss gainLoss : gainsLosses) {
            if (gainLoss.chromosome().equals(chromosome) && gainLoss.chromosomeBand().equals(chromosomeBand)) {
                return true;
            }
        }

        return false;
    }
}
