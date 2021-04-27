package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.drivercatalog.CNADrivers;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.BestFit;
import com.hartwig.hmftools.common.purple.PurpleQC;

import org.jetbrains.annotations.NotNull;

public final class PurpleQCFactory
{
    private PurpleQCFactory()
    {
    }

    @NotNull
    public static PurpleQC create(double contamination, @NotNull BestFit bestFit, @NotNull Gender amberGender, @NotNull Gender cobaltGender,
            @NotNull List<PurpleCopyNumber> copyNumbers, @NotNull List<GeneCopyNumber> geneCopyNumbers,
            @NotNull final Set<GermlineAberration> aberrations, int amberMeanDepth)
    {
        boolean containsAnySvSupport = copyNumbers.stream().anyMatch(PurpleCopyNumber::svSupport);

        int unsupportedCopyNumberSegments = containsAnySvSupport ? (int) copyNumbers.stream()
                .filter(x -> x.segmentStartSupport() == NONE && x.segmentEndSupport() == NONE)
                .count() : 0;
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
                .amberMeanDepth(amberMeanDepth)
                .build();
    }
}
