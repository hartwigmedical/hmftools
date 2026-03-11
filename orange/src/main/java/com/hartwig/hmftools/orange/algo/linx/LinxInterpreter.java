package com.hartwig.hmftools.orange.algo.linx;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Objects;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.PurpleAllelicDepth;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;
import com.hartwig.hmftools.orange.conversion.ConversionUtil;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.jetbrains.annotations.Nullable;

public class LinxInterpreter
{
    private final CytoBands mCytoBands;

    public LinxInterpreter(final CytoBands cytoBands)
    {
        mCytoBands = cytoBands;
    }

    public LinxRecord interpret(final LinxData linx, @Nullable final IsofoxData isofoxData)
    {
        LOGGER.info("Analysing Linx data");

        LinxBreakendInterpreter somaticBreakendInterpreter = new LinxBreakendInterpreter(
                linx.somaticSvAnnotations(), linx.somaticDrivers(), mCytoBands);

        LinxBreakendInterpreter germlineBreakendInterpreter = new LinxBreakendInterpreter(
                Objects.requireNonNullElse(linx.germlineSvAnnotations(), List.of()), linx.germlineDrivers(), mCytoBands);

        return ImmutableLinxRecord.builder()
                .somaticStructuralVariants(ConversionUtil.mapToIterable(linx.somaticSvAnnotations(), LinxConversion::convert))
                .somaticDrivers(ConversionUtil.mapToIterable(linx.somaticDriverData(), LinxConversion::convert))
                .fusions(buildFusions(linx.fusions(), linx.somaticBreakends(), isofoxData))
                .somaticBreakends(somaticBreakendInterpreter.convertBreakends(linx.somaticBreakends()))
                .germlineStructuralVariants(ConversionUtil.mapToIterable(linx.germlineSvAnnotations(), LinxConversion::convert))
                .germlineBreakends(germlineBreakendInterpreter.convertBreakends(linx.germlineBreakends()))
                .build();
    }

    private static List<LinxFusion> buildFusions(
            final List<com.hartwig.hmftools.common.linx.LinxFusion> fusions, final List<LinxBreakend> breakends,
            @Nullable final IsofoxData isofoxData)
    {
        List<LinxFusion> convertedFusions = Lists.newArrayListWithCapacity(fusions.size());

        for(com.hartwig.hmftools.common.linx.LinxFusion fusion : fusions)
        {
            LinxBreakend breakendUp = breakends.stream().filter(x -> x.id() == fusion.fivePrimeBreakendId()).findFirst().orElse(null);
            LinxBreakend breakendDown = breakends.stream().filter(x -> x.id() == fusion.threePrimeBreakendId()).findFirst().orElse(null);

            if(breakendUp == null || breakendDown == null)
            {
                LOGGER.error("fusion({}) missing corresponding breakends", fusion);
                continue;
            }

            LinxFusionType fusionType = LinxFusionType.valueOf(fusion.reportedType().toString());

            String contextUp = buildContextStr(breakendUp.regionType(), fusionType, fusion.fusedExonUp());
            String contextDown = buildContextStr(breakendDown.regionType(), fusionType, fusion.fusedExonDown());

            double avgJcn = (breakendUp.undisruptedCopyNumber() + breakendDown.undisruptedCopyNumber()) * 0.5;

            PurpleAllelicDepth rnaSupport = findRnaSupport(fusion, isofoxData);

            LinxFusion convertedFusion = ImmutableLinxFusion.builder()
                    .geneUp(breakendUp.gene())
                    .contextUp(contextUp)
                    .transcriptUp(breakendUp.transcriptId())
                    .geneDown(breakendDown.gene())
                    .contextDown(contextDown)
                    .transcriptDown(breakendDown.transcriptId())
                    .reportedType(LinxFusionType.valueOf(fusion.reportedType()))
                    .phased(FusionPhasedType.valueOf(fusion.phased().name()))
                    .driverInterpretation(fromFusionLikelihood(fusion.likelihood()))
                    .fusedExonUp(fusion.fusedExonUp())
                    .fusedExonDown(fusion.fusedExonDown())
                    .chainLinks(fusion.chainLinks())
                    .chainTerminated(fusion.chainTerminated())
                    .domainsKept(fusion.domainsKept())
                    .domainsLost(fusion.domainsLost())
                    .junctionCopyNumber(avgJcn)
                    .rnaSupport(rnaSupport)
                    .build();


            convertedFusions.add(convertedFusion);
        }

        return convertedFusions;
    }

    private static PurpleAllelicDepth findRnaSupport(com.hartwig.hmftools.common.linx.LinxFusion fusion, final IsofoxData isofoxData)
    {
        if(isofoxData == null)
            return null;

        int totalCoverage = 0;
        int totalFragments = 0;

        for(RnaFusion rnaFusion : isofoxData.fusions())
        {
            if(rnaFusion.name().equals(fusion.name()))
            {
                int fragments = rnaFusion.splitFragments();
                int averageDepth = min((int)round((rnaFusion.depthUp() + rnaFusion.depthDown()) * 0.5), fragments);
                totalCoverage += averageDepth;
                totalFragments += fragments;
            }
        }

        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(totalFragments)
                .totalReadCount(totalCoverage)
                .build();
    }

    private static DriverInterpretation fromFusionLikelihood(final FusionLikelihoodType fusionLikelihoodType)
    {
        switch(fusionLikelihoodType)
        {
            case HIGH: return DriverInterpretation.HIGH;
            case LOW: return DriverInterpretation.LOW;
            default: return DriverInterpretation.UNKNOWN;
        }
    }

    private static String buildContextStr(final TranscriptRegionType regionType, final LinxFusionType knownType, int fusedExon)
    {
        switch(regionType)
        {
            case UPSTREAM:
                return "Promoter Region";
            case DOWNSTREAM:
                return "Post-coding";
            case IG:
                return "IG";
            case EXONIC:
            case INTRONIC:
                return String.format("Exon %d", fusedExon);
            case UNKNOWN:
            {
                if(knownType == LinxFusionType.PROMISCUOUS_ENHANCER_TARGET)
                {
                    return "Unknown";
                }
                else
                {
                    return String.format("ERROR: %s", regionType);
                }
            }
        }

        throw new IllegalStateException("TranscriptRegionType not supported in determination of fusion context: " + regionType);
    }
}
