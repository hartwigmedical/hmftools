package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;
import static com.hartwig.hmftools.orange.algo.util.DriverUtils.convertReportedStatus;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;

import org.jetbrains.annotations.Nullable;

public class LinxBreakendInterpreter
{
    private final Map<Integer, LinxSvAnnotation> mLinxSvAnnotationsMap;
    private final EnsemblDataCache mEnsemblDataCache;

    public LinxBreakendInterpreter(final List<LinxSvAnnotation> linxSvAnnotations, final EnsemblDataCache ensemblDataCache)
    {
        mLinxSvAnnotationsMap = linxSvAnnotations.stream().collect(Collectors.toMap(LinxSvAnnotation::svId, s -> s));
        mEnsemblDataCache = ensemblDataCache;
    }

    public LinxBreakend interpret(com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        LinxSvAnnotation svAnnotation = mLinxSvAnnotationsMap.get(linxBreakend.svId());
        if(svAnnotation == null)
        {
            // CHECK: can't see why this would ever occur or be acceptable it if did
            LOGGER.warn("No Linx SV annotation found for breakend {}", linxBreakend.id());
        }

        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .chromosome(chromosome(svAnnotation, linxBreakend.isStart()))
                .chromosomeBand(chromosomeBand(linxBreakend.gene()))
                .transcript(linxBreakend.transcriptId())
                .isCanonical(linxBreakend.canonical())
                .geneOrientation(geneOrientation(linxBreakend.geneOrientation()))
                .isCanonical(linxBreakend.canonical())
                .disruptive(linxBreakend.disruptive())
                .reportedStatus(convertReportedStatus(linxBreakend.reportedStatus()))
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .type(svAnnotation != null ? LinxBreakendType.valueOf(svAnnotation.type().name()) : LinxBreakendType.BND)
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .orientation(orientation(svAnnotation, linxBreakend.isStart()))
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                .junctionCopyNumber(junctionCopyNumber(svAnnotation))
                .build();
    }

    private static String chromosome(final LinxSvAnnotation svAnnotation, boolean isStart)
    {
        if(svAnnotation == null)
            return "";

        String coords = isStart ? svAnnotation.coordsStart() : svAnnotation.coordsEnd();
        return com.hartwig.hmftools.common.linx.LinxBreakend.chromosomeFromCoords(coords);
    }

    private static byte orientation(final LinxSvAnnotation svAnnotation, boolean isStart)
    {
        if(svAnnotation == null)
            return 0;

        String coords = isStart ? svAnnotation.coordsStart() : svAnnotation.coordsEnd();
        return com.hartwig.hmftools.common.linx.LinxBreakend.orientationFromCoords(coords);
    }

    private String chromosomeBand(final String gene)
    {
        GeneData geneData = mEnsemblDataCache.getGeneDataByName(gene);
        return geneData != null ? geneData.KaryotypeBand : "";
    }

    private LinxGeneOrientation geneOrientation(final String geneOrientation)
    {
        if(geneOrientation.equals(com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM))
        {
            return LinxGeneOrientation.UPSTREAM;
        }
        else if(geneOrientation.equals(com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_DOWNSTREAM))
        {
            return LinxGeneOrientation.DOWNSTREAM;
        }
        throw new IllegalArgumentException("Unknown gene orientation: " + geneOrientation);
    }

    @VisibleForTesting
    static double junctionCopyNumber(@Nullable LinxSvAnnotation svAnnotation)
    {
        if(svAnnotation == null)
        {
            return 0;
        }
        else
        {
            return 0.5D * (svAnnotation.junctionCopyNumberMin() + svAnnotation.junctionCopyNumberMax());
        }
    }
}