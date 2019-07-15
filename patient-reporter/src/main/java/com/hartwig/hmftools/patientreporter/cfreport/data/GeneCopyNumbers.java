package com.hartwig.hmftools.patientreporter.cfreport.data;

import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberInterpretation;

import org.jetbrains.annotations.NotNull;

public final class GeneCopyNumbers {

    private GeneCopyNumbers() {
    }

    @NotNull
    public static List<GeneCopyNumber> sort(@NotNull final List<GeneCopyNumber> reportableCopyNumbers) {
        return reportableCopyNumbers.stream().sorted((copyNumber1, copyNumber2) -> {
            String location1 = GeneUtil.zeroPrefixed(copyNumber1.chromosome() + copyNumber1.chromosomeBand());
            String location2 = GeneUtil.zeroPrefixed(copyNumber2.chromosome() + copyNumber2.chromosomeBand());

            if (location1.equals(location2)) {
                return copyNumber1.gene().compareTo(copyNumber2.gene());
            } else {
                return location1.compareTo(location2);
            }
        }).collect(Collectors.toList());
    }

    @NotNull
    public static Set<String> amplifiedGenes(@NotNull List<GeneCopyNumber> reportableCopyNumbers) {
        final Set<String> genes = Sets.newHashSet();
        for (GeneCopyNumber copyNumber : reportableCopyNumbers) {
            if (interpret(copyNumber) == CopyNumberInterpretation.GAIN) {
                genes.add(copyNumber.gene());
            }
        }
        return genes;
    }

    @NotNull
    public static Set<String> lossGenes(@NotNull List<GeneCopyNumber> reportableCopyNumbers) {
        final Set<String> genes = Sets.newHashSet();
        for (GeneCopyNumber copyNumber : reportableCopyNumbers) {
            if (interpret(copyNumber) == CopyNumberInterpretation.LOSS) {
                genes.add(copyNumber.gene());
            }
        }
        return genes;
    }

    @NotNull
    public static String type(@NotNull GeneCopyNumber geneCopyNumber) {
        CopyNumberInterpretation interpretation = interpret(geneCopyNumber);
        if (interpretation == CopyNumberInterpretation.GAIN) {
            return "gain";
        } else {
            // At this point we only have losses and gains.
            assert interpretation == CopyNumberInterpretation.LOSS;
            if (geneCopyNumber.maxCopyNumber() < 0.5) {
                return "full loss";
            } else {
                return "partial loss";
            }
        }
    }

    @NotNull
    private static CopyNumberInterpretation interpret(@NotNull GeneCopyNumber geneCopyNumber) {
        return CopyNumberInterpretation.fromCopyNumber(geneCopyNumber.minCopyNumber());
    }
}
