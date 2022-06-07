package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.ReportableGainLoss;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class CopyNumberSelector {

    private CopyNumberSelector() {
    }

    @NotNull
    public static List<ReportableGainLoss> selectNonDriverGains(@NotNull List<ReportableGainLoss> unreportedGainsLosses) {
        List<ReportableGainLoss> unreportedFullGains = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.FULL_GAIN)
                .collect(Collectors.toList());

        Map<CopyNumberKey, ReportableGainLoss> maxGainPerLocation = Maps.newHashMap();
        for (ReportableGainLoss gain : unreportedFullGains) {
            CopyNumberKey key = new CopyNumberKey(gain.chromosome(), gain.chromosomeBand());
            ReportableGainLoss maxGain = maxGainPerLocation.get(key);
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
    public static List<ReportableGainLoss> selectNonDriverLosses(@NotNull List<ReportableGainLoss> unreportedGainsLosses,
            @NotNull List<ReportableGainLoss> reportedGainLosses) {
        List<ReportableGainLoss> lossesNoAllosomes = Lists.newArrayList();

        List<ReportableGainLoss> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        List<ReportableGainLoss> reportedLosses = reportedGainLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        for (ReportableGainLoss loss : unreportedLosses) {
            if (!HumanChromosome.fromString(loss.chromosome()).isAllosome() && !locusPresent(reportedLosses,
                    loss.chromosome(),
                    loss.chromosomeBand())) {
                lossesNoAllosomes.add(loss);
            }
        }

        Map<CopyNumberKey, ReportableGainLoss> oneLossPerLocation = Maps.newHashMap();
        for (ReportableGainLoss loss : lossesNoAllosomes) {
            CopyNumberKey key = new CopyNumberKey(loss.chromosome(), loss.chromosomeBand());
            if (!oneLossPerLocation.containsKey(key)) {
                oneLossPerLocation.put(key, loss);
            }
        }

        return Lists.newArrayList(oneLossPerLocation.values().iterator());
    }

    private static boolean locusPresent(@NotNull List<ReportableGainLoss> gainsLosses, @NotNull String chromosome,
            @NotNull String chromosomeBand) {
        for (ReportableGainLoss gainLoss : gainsLosses) {
            if (gainLoss.chromosome().equals(chromosome) && gainLoss.chromosomeBand().equals(chromosomeBand)) {
                return true;
            }
        }

        return false;
    }
}
