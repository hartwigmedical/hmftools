package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;

public class LinxBreakendInterpreter
{
    private final Map<Integer, LinxSvAnnotation> mLinxSvAnnotationsMap;
    private final CytoBands mCytoBands;

    public LinxBreakendInterpreter(final List<LinxSvAnnotation> linxSvAnnotations, final CytoBands cytoBands)
    {
        mLinxSvAnnotationsMap = linxSvAnnotations.stream().collect(Collectors.toMap(LinxSvAnnotation::svId, s -> s));
        mCytoBands = cytoBands;
    }

    public LinxBreakend interpret(com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        LinxSvAnnotation svAnnotation = mLinxSvAnnotationsMap.get(linxBreakend.svId());

        String breakendCoords = linxBreakend.isStart() ? svAnnotation.coordsStart() : svAnnotation.coordsEnd();

        String chromosome = com.hartwig.hmftools.common.linx.LinxBreakend.chromosomeFromCoords(breakendCoords);
        int position = com.hartwig.hmftools.common.linx.LinxBreakend.positionFromCoords(breakendCoords);
        byte orientation = com.hartwig.hmftools.common.linx.LinxBreakend.orientationFromCoords(breakendCoords);

        String cytoBand = mCytoBands.getCytoBandName(chromosome, position);

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
                .build();
    }

    @VisibleForTesting
    static double junctionCopyNumber(final LinxSvAnnotation svAnnotation)
    {
        return 0.5D * (svAnnotation.junctionCopyNumberMin() + svAnnotation.junctionCopyNumberMax());
    }
}