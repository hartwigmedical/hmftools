package com.hartwig.hmftools.orange.algo.linx;

import static java.lang.Math.min;
import static java.lang.Math.round;

import static com.hartwig.hmftools.common.linx.LinxCommonTypes.formVisPlotClusterPrefix;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.LinxBreakend;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.rna.RnaFusion;
import com.hartwig.hmftools.datamodel.common.AllelicDepth;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.common.ImmutableAllelicDepth;
import com.hartwig.hmftools.orange.algo.isofox.IsofoxData;

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
        LOGGER.debug("analysing Linx data");

        List<com.hartwig.hmftools.datamodel.linx.LinxBreakend> somaticBreakends = LinxBreakendInterpreter.buildSomaticBreakends(
                linx, mCytoBands);

        List<com.hartwig.hmftools.datamodel.linx.LinxBreakend> germlineBreakends = LinxBreakendInterpreter.buildGermlineBreakends(
                linx, mCytoBands);

        return ImmutableLinxRecord.builder()
                .fusions(buildFusions(linx, isofoxData))
                .somaticBreakends(somaticBreakends)
                .germlineBreakends(germlineBreakends)
                .build();
    }

    private static LinxBreakend findBreakend(final int breakendId, final List<LinxBreakend> breakends)
    {
        return breakends.stream().filter(x -> x.id() == breakendId).findFirst().orElse(null);
    }

    private static LinxSvAnnotation findSvAnnotation(final LinxBreakend breakend, final List<LinxSvAnnotation> svAnnotations)
    {
        return svAnnotations.stream().filter(x -> x.svId() == breakend.svId()).findFirst().orElse(null);
    }

    public static String findReportableLinxPlot(final List<String> linxPlots, final int clusterId)
    {
        String plotPrefix = formVisPlotClusterPrefix(String.valueOf(clusterId));
        return linxPlots.stream().filter(x -> x.contains(plotPrefix)).findFirst().orElse(null);
    }

    private static List<LinxFusion> buildFusions(final LinxData linx, @Nullable final IsofoxData isofoxData)
    {
        List<LinxFusion> convertedFusions = Lists.newArrayListWithCapacity(linx.fusions().size());

        for(com.hartwig.hmftools.common.linx.LinxFusion fusion : linx.fusions())
        {
            LinxBreakend breakendUp = findBreakend(fusion.fivePrimeBreakendId(), linx.somaticBreakends());
            LinxBreakend breakendDown = findBreakend(fusion.threePrimeBreakendId(), linx.somaticBreakends());

            if(breakendUp == null || breakendDown == null)
            {
                LOGGER.error("fusion({}) missing corresponding breakends", fusion);
                continue;
            }

            LinxSvAnnotation svAnnotation = findSvAnnotation(breakendUp, linx.somaticSvAnnotations());

            if(svAnnotation == null)
            {
                LOGGER.error("fusion({}) missing corresponding svAnnotation", fusion);
                continue;
            }

            String plotFilename = findReportableLinxPlot(linx.reportableEventPlots(), svAnnotation.clusterId());

            LinxFusionType fusionType = LinxFusionType.valueOf(fusion.reportedType().toString());

            String geneUp = breakendUp.gene();
            String geneDown = breakendDown.gene();

            if(fusionType == LinxFusionType.ENHANCER_KNOWN_PAIR || fusionType == LinxFusionType.ENHANCER_PROMISCUOUS)
            {
                geneUp = breakendUp.transcriptId(); // to use the IG / TC name
            }

            String contextUp = buildContextStr(breakendUp.regionType(), fusionType, fusion.fusedExonUp());
            String contextDown = buildContextStr(breakendDown.regionType(), fusionType, fusion.fusedExonDown());

            double avgJcn = (breakendUp.undisruptedCopyNumber() + breakendDown.undisruptedCopyNumber()) * 0.5;

            AllelicDepth rnaSupport = findRnaSupport(fusion, isofoxData);

            LinxFusion convertedFusion = ImmutableLinxFusion.builder()
                    .geneUp(geneUp)
                    .contextUp(contextUp)
                    .transcriptUp(breakendUp.transcriptId())
                    .geneDown(geneDown)
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
                    .plotFilename(plotFilename)
                    .build();

            convertedFusions.add(convertedFusion);
        }

        return convertedFusions;
    }

    private static AllelicDepth findRnaSupport(com.hartwig.hmftools.common.linx.LinxFusion fusion, final IsofoxData isofoxData)
    {
        if(isofoxData == null)
            return null;

        int totalCoverage = 0;
        int totalFragments = 0;
        boolean matched = false;

        for(RnaFusion rnaFusion : isofoxData.fusions())
        {
            if(rnaFusion.name().equals(fusion.name()))
            {
                matched = true;

                int fragments = rnaFusion.splitFragments() + rnaFusion.discordantFrags() + rnaFusion.realignedFrags();
                int averageDepth = min((int)round((rnaFusion.depthUp() + rnaFusion.depthDown()) * 0.5), fragments);
                totalCoverage += averageDepth;
                totalFragments += fragments;
            }
        }

        if(!matched)
            return null;

        return ImmutableAllelicDepth.builder()
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
