package com.hartwig.hmftools.purple;

import static com.hartwig.hmftools.common.purple.SegmentSupport.NONE;
import static com.hartwig.hmftools.purple.fitting.FittedPurityScoreFactory.polyclonalProportion;
import static com.hartwig.hmftools.purple.fitting.WholeGenomeDuplication.wholeGenomeDuplication;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.purple.drivers.DeletionDrivers;
import com.hartwig.hmftools.common.genome.chromosome.GermlineAberration;
import com.hartwig.hmftools.common.purple.FittedPurityMethod;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.ImmutablePurityContext;
import com.hartwig.hmftools.common.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
import com.hartwig.hmftools.common.purple.TumorMutationalStatus;
import com.hartwig.hmftools.common.purple.MicrosatelliteStatus;
import com.hartwig.hmftools.purple.copynumber.LohCalcData;
import com.hartwig.hmftools.purple.copynumber.LohCalcs;
import com.hartwig.hmftools.purple.fitting.BestFit;
import com.hartwig.hmftools.purple.somatic.SomaticStream;
import com.hartwig.hmftools.purple.sv.SomaticSvCache;

public final class PurpleSummaryData
{
    public static PurpleQC createQC(
            double contamination, final BestFit bestFit, final Gender amberGender, final Gender cobaltGender,
            final List<PurpleCopyNumber> copyNumbers, final List<GeneCopyNumber> geneCopyNumbers,
            final Set<GermlineAberration> aberrations, int amberMeanDepth, int maxDeletedGenes, double tincLevel, final double chimerismPercentage)
    {
        boolean containsAnySvSupport = copyNumbers.stream().anyMatch(PurpleCopyNumber::svSupport);

        int unsupportedCopyNumberSegments = containsAnySvSupport ? (int) copyNumbers.stream()
                .filter(x -> x.segmentStartSupport() == NONE && x.segmentEndSupport() == NONE)
                .count() : 0;

        int deletedGenes = DeletionDrivers.deletedGenes(geneCopyNumbers);

        Set<PurpleQCStatus> statusSet = PurpleQCStatus.calcStatus(
                PurpleQCStatus.genderPass(amberGender, cobaltGender, aberrations),
                unsupportedCopyNumberSegments, deletedGenes, bestFit.Fit.purity(), bestFit.Method, contamination, maxDeletedGenes, tincLevel);

        LohCalcData lohCalcData = LohCalcs.calcLohData(copyNumbers);

        return ImmutablePurpleQC.builder()
                .status(statusSet)
                .method(bestFit.Method)
                .contamination(contamination)
                .cobaltGender(cobaltGender)
                .amberGender(amberGender)
                .purity(bestFit.Fit.purity())
                .copyNumberSegments(copyNumbers.size())
                .unsupportedCopyNumberSegments(unsupportedCopyNumberSegments)
                .deletedGenes(deletedGenes)
                .germlineAberrations(aberrations)
                .amberMeanDepth(amberMeanDepth)
                .lohPercent(lohCalcData.lohPercent())
                .tincLevel(tincLevel)
                .chimerismPercentage(chimerismPercentage)
                .build();
    }

    public static PurityContext createPurity(
            final BestFit bestFit, final Gender gender, final PurpleConfig config, final PurpleQC qcChecks,
            final List<PurpleCopyNumber> copyNumbers, final SomaticStream somaticStream, final SomaticSvCache svCache)
    {
        boolean wholeGenomeDuplication = !config.tumorOnlyMode() || bestFit.Method == FittedPurityMethod.NORMAL
                ? wholeGenomeDuplication(copyNumbers) : false;

        return ImmutablePurityContext.builder()
                .bestFit(bestFit.Fit)
                .method(bestFit.Method)
                .gender(gender)
                .runMode(config.runMode())
                .score(bestFit.Score)
                .polyClonalProportion(polyclonalProportion(copyNumbers))
                .wholeGenomeDuplication(wholeGenomeDuplication)
                .microsatelliteIndelsPerMb(somaticStream != null ? somaticStream.msiIndelsPerMb() : 0)
                .microsatelliteStatus(somaticStream != null ? somaticStream.microsatelliteStatus() : MicrosatelliteStatus.UNKNOWN)
                .tumorMutationalLoad(somaticStream != null ? somaticStream.tumorMutationalLoad() : 0)
                .tumorMutationalLoadStatus(
                        somaticStream != null ? somaticStream.tumorMutationalLoadStatus() : TumorMutationalStatus.UNKNOWN)
                .tumorMutationalBurdenPerMb(somaticStream != null ? somaticStream.tumorMutationalBurdenPerMb() : 0)
                .tumorMutationalBurdenStatus(
                        somaticStream != null ? somaticStream.tumorMutationalBurdenPerMbStatus() : TumorMutationalStatus.UNKNOWN)
                .svTumorMutationalBurden(svCache.passingBnd())
                .qc(qcChecks)
                .targeted(config.TargetRegionsMode)
                .build();
    }
}
