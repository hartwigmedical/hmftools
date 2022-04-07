package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.segment.SegmentSupport.NONE;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.wholeGenomeDuplication;
import static com.hartwig.hmftools.purple.purity.FittedPurityScoreFactory.polyclonalProportion;

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
import com.hartwig.hmftools.common.purple.purity.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus;
import com.hartwig.hmftools.common.variant.tml.TumorMutationalStatus;
import com.hartwig.hmftools.purple.config.PurpleConfig;
import com.hartwig.hmftools.purple.somatic.SomaticStream;

public final class PurpleSummaryData
{
    public static PurpleQC createQC(double contamination, final BestFit bestFit, final Gender amberGender, final Gender cobaltGender,
            final List<PurpleCopyNumber> copyNumbers, final List<GeneCopyNumber> geneCopyNumbers,
            final Set<GermlineAberration> aberrations, int amberMeanDepth)
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

    public static PurityContext createPurity(
            final String version, final BestFit bestFit, final Gender gender, final PurpleConfig config, final PurpleQC qcChecks,
            final List<PurpleCopyNumber> copyNumbers, final SomaticStream somaticStream, final StructuralVariantCache svCache)
    {
        return ImmutablePurityContext.builder()
                .version(version)
                .bestFit(bestFit.fit())
                .method(bestFit.method())
                .gender(gender)
                .runMode(config.runMode())
                .score(bestFit.score())
                .polyClonalProportion(polyclonalProportion(copyNumbers))
                .wholeGenomeDuplication(wholeGenomeDuplication(copyNumbers))
                .microsatelliteIndelsPerMb(somaticStream != null ? somaticStream.msiIndelsPerMb() : 0)
                .microsatelliteStatus(somaticStream != null ? somaticStream.microsatelliteStatus() : MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalLoad(somaticStream != null ? somaticStream.tumorMutationalLoad() : 0)
                .tumorMutationalLoadStatus(somaticStream != null ? somaticStream.tumorMutationalLoadStatus() : TumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(somaticStream != null ? somaticStream.tumorMutationalBurdenPerMb() : 0)
                .tumorMutationalBurdenStatus(somaticStream != null ? somaticStream.tumorMutationalBurdenPerMbStatus() : TumorMutationalStatus.UNKNOWN)
                .svTumorMutationalBurden(svCache.passingBnd())
                .qc(qcChecks)
                .targeted(config.targetRegionsMode())
                .build();
    }
}
