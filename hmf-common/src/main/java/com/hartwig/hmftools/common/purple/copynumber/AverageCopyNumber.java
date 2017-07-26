package com.hartwig.hmftools.common.purple.copynumber;

import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.region.GenomeRegion;

import org.jetbrains.annotations.NotNull;

public class AverageCopyNumber {

    public static double averageCopyNumber(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        final String key = "key";
        final List<PurpleCopyNumber> filtered =
                copyNumbers.stream().filter(x -> !x.chromosome().equals("X") && !x.chromosome().equals("Y")).collect(Collectors.toList());

        return averageCopyNumber(x -> key, filtered).getOrDefault(key, 0D);
    }

    public static Map<String, Double> averageChromosomalCopyNumber(@NotNull final List<PurpleCopyNumber> copyNumbers) {
        return averageCopyNumber(GenomeRegion::chromosome, copyNumbers);
    }

    private static Map<String, Double> averageCopyNumber(Function<PurpleCopyNumber, String> keyFunction,
            @NotNull final List<PurpleCopyNumber> copyNumbers) {

        final Map<String, Double> result = Maps.newHashMap();
        final Map<String, Double> sumMap = Maps.newHashMap();
        final Map<String, Integer> bafCountMap = Maps.newHashMap();

        for (PurpleCopyNumber copyNumber : copyNumbers) {
            final String key = keyFunction.apply(copyNumber);
            final int bafCount = copyNumber.bafCount();
            final double weightedCopyNumber = copyNumber.averageTumorCopyNumber() * bafCount;
            bafCountMap.merge(key, bafCount, (aDouble, aDouble2) -> aDouble + aDouble2);
            sumMap.merge(key, weightedCopyNumber, (aDouble, aDouble2) -> aDouble + aDouble2);
        }

        for (String key : bafCountMap.keySet()) {
            int bafCount = bafCountMap.getOrDefault(key, 0);
            double copyNumberSum = sumMap.getOrDefault(key, 0D);
            double average = bafCount == 0 ? 0 : copyNumberSum / bafCount;
            result.put(key, average);
        }

        return result;
    }
}
