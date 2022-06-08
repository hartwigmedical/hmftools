package com.hartwig.hmftools.orange.algo.selection;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.common.purple.interpretation.CopyNumberInterpretation;
import com.hartwig.hmftools.common.purple.interpretation.GainLoss;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.NotNull;

public final class CopyNumberSelector {

    private CopyNumberSelector() {
    }

    @NotNull
    public static List<GainLoss> selectNonDriverGains(@NotNull List<GainLoss> unreportedGainsLosses) {
        List<GainLoss> unreportedFullGains = unreportedGainsLosses.stream()
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
    public static List<GainLoss> selectNonDriverLosses(@NotNull List<GainLoss> unreportedGainsLosses,
            @NotNull List<GainLoss> reportedGainLosses) {
        List<GainLoss> lossesNoAllosomes = Lists.newArrayList();

        List<GainLoss> unreportedLosses = unreportedGainsLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        List<GainLoss> reportedLosses = reportedGainLosses.stream()
                .filter(gainLoss -> gainLoss.interpretation() == CopyNumberInterpretation.PARTIAL_LOSS
                        || gainLoss.interpretation() == CopyNumberInterpretation.FULL_LOSS)
                .collect(Collectors.toList());

        for (GainLoss loss : unreportedLosses) {
            if (!HumanChromosome.fromString(loss.chromosome()).isAllosome() && !locusPresent(reportedLosses,
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

    private static boolean locusPresent(@NotNull List<GainLoss> gainsLosses, @NotNull String chromosome,
            @NotNull String chromosomeBand) {
        for (GainLoss gainLoss : gainsLosses) {
            if (gainLoss.chromosome().equals(chromosome) && gainLoss.chromosomeBand().equals(chromosomeBand)) {
                return true;
            }
        }

        return false;
    }
}
