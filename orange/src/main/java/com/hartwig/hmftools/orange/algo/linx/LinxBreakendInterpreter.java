package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.algo.linx.LinxInterpreter.findReportableLinxPlot;

import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverType;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxDriverType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;

public final class LinxBreakendInterpreter
{
    public static List<LinxBreakend> buildSomaticBreakends(final LinxData linxData, final CytoBands cytoBands)
    {
        List<LinxBreakend> reportedBreakends = Lists.newArrayList();

        Map<Integer,LinxSvAnnotation> svAnnotationsMap = linxData.somaticSvAnnotations().stream()
                .collect(Collectors.toMap(LinxSvAnnotation::svId, s -> s));

        for(com.hartwig.hmftools.common.linx.LinxBreakend breakend : linxData.somaticBreakends())
        {
            if(breakend.reportedStatus() != com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED)
                continue;

            LinxBreakend convertedBreakend = build(
                    breakend, svAnnotationsMap, linxData.somaticDrivers(), cytoBands, linxData.reportableEventPlots());

            reportedBreakends.add(convertedBreakend);
        }

        return reportedBreakends;
    }

    public static List<LinxBreakend> buildGermlineBreakends(final LinxData linxData, final CytoBands cytoBands)
    {
        List<LinxBreakend> reportedBreakends = Lists.newArrayList();

        Map<Integer,LinxSvAnnotation> svAnnotationsMap = linxData.germlineSvAnnotations().stream()
                .collect(Collectors.toMap(LinxSvAnnotation::svId, s -> s));

        for(com.hartwig.hmftools.common.linx.LinxBreakend breakend : linxData.germlineBreakends())
        {
            if(breakend.reportedStatus() != com.hartwig.hmftools.common.purple.ReportedStatus.REPORTED)
                continue;

            LinxBreakend convertedBreakend = build(
                    breakend, svAnnotationsMap, linxData.germlineDrivers(), cytoBands, Collections.emptyList());

            reportedBreakends.add(convertedBreakend);
        }

        return reportedBreakends;
    }

    public static LinxBreakend build(
            final com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend, final Map<Integer, LinxSvAnnotation> svAnnotationsMap,
            final List<DriverCatalog> drivers, final CytoBands cytoBands, final List<String> linxPlots)
    {
        LinxSvAnnotation svAnnotation = svAnnotationsMap.get(linxBreakend.svId());

        String breakendCoords = linxBreakend.isStart() ? svAnnotation.coordsStart() : svAnnotation.coordsEnd();

        String chromosome = com.hartwig.hmftools.common.linx.LinxBreakend.chromosomeFromCoords(breakendCoords);
        int position = com.hartwig.hmftools.common.linx.LinxBreakend.positionFromCoords(breakendCoords);
        byte orientation = com.hartwig.hmftools.common.linx.LinxBreakend.orientationFromCoords(breakendCoords);

        String cytoBand = cytoBands.getCytoBandName(chromosome, position);

        DriverCatalog driverCatalog = drivers.stream().filter(x -> x.gene().equals(linxBreakend.gene())).findFirst().orElse(null);

        LinxDriverType driverType = LinxDriverType.DISRUPTION;
        double driverLikelihood = 0;

        if(driverCatalog != null)
        {
            if(driverCatalog.driver() != DriverType.GERMLINE_DISRUPTION)
                driverType = LinxDriverType.valueOf(driverCatalog.driver().toString());

            driverLikelihood = driverCatalog.driverLikelihood();
        }

        String plotFilename = findReportableLinxPlot(linxPlots, svAnnotation.clusterId());

        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .chromosome(chromosome)
                .chromosomeBand(cytoBand)
                .transcript(linxBreakend.transcriptId())
                .isCanonical(linxBreakend.canonical())
                .geneOrientation(linxBreakend.geneOrientation().equals(com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM) ?
                        LinxGeneOrientation.UPSTREAM : LinxGeneOrientation.DOWNSTREAM)
                .disruptive(linxBreakend.disruptive())
                .reportedStatus(ReportedStatus.valueOf(linxBreakend.reportedStatus().name()))
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .type(svAnnotation != null ? LinxBreakendType.valueOf(svAnnotation.type().name()) : LinxBreakendType.BND)
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .orientation(orientation)
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                .junctionCopyNumber(junctionCopyNumber(svAnnotation))
                .driverType(driverType)
                .driverLikelihood(driverLikelihood)
                .plotFilename(plotFilename)
                .build();
    }

    @VisibleForTesting
    static double junctionCopyNumber(final LinxSvAnnotation svAnnotation)
    {
        return 0.5D * (svAnnotation.junctionCopyNumberMin() + svAnnotation.junctionCopyNumberMax());
    }
}