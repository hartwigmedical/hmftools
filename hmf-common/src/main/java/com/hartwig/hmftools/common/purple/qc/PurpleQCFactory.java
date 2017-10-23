package com.hartwig.hmftools.common.purple.qc;

import java.util.List;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.segment.StructuralVariantSupport;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFactory {

    @NotNull
    public static PurpleQC create(@NotNull FittedPurity purity, @NotNull List<PurpleCopyNumber> copyNumbers, @NotNull Gender purpleGender, @NotNull Gender cobaltGender) {
        final List<PurpleCopyNumber> trailingSegments = copyNumbers.stream().filter(x -> x.start() != 1L).collect(Collectors.toList());
        int ratioOnlySegments = (int) trailingSegments.stream().filter(x -> x.structuralVariantSupport() == StructuralVariantSupport.NONE).count();

        return ImmutablePurpleQC.builder()
                .cobaltGender(cobaltGender)
                .amberGender(purpleGender)
                .ploidy(purity.ploidy())
                .ratioSegments(ratioOnlySegments)
                .trailingSegments(trailingSegments.size())
                .build();
    }
}
