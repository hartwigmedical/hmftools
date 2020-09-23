package com.hartwig.hmftools.common.purple.qc;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.CENTROMERE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.INF;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.TELOMERE;

import java.util.EnumSet;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.segment.SegmentSupport;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFactory {

    private static final EnumSet<SegmentSupport> INTERNAL = EnumSet.of(NONE, CENTROMERE, TELOMERE, INF);

    private PurpleQCFactory() {
    }

    @NotNull
    public static PurpleQC create(double contamination, @NotNull BestFit bestFit, @NotNull Gender amberGender, @NotNull Gender cobaltGender,
            @NotNull List<PurpleCopyNumber> copyNumbers, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull final Set<GermlineAberration> aberrations) {

        int unsupportedCopyNumberSegments = (int) copyNumbers.stream().filter(x -> !x.svSupport()).count();
        int deletedGenes = CNADrivers.deletedGenes(geneCopyNumbers);

        return ImmutablePurpleQC.builder()
                .method(bestFit.method())
                .contamination(contamination)
                .cobaltGender(cobaltGender)
                .amberGender(amberGender)
                .purity(bestFit.fit().purity())
                .copyNumberSegments(copyNumbers.size())
                .unsupportedCopyNumberSegments(unsupportedCopyNumberSegments)
                .deletedGenes(deletedGenes)
                .germlineAberrations(aberrations)
                .build();
    }
}
