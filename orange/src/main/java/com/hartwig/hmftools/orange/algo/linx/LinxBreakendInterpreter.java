package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.common.linx.LinxBreakend.BREAKEND_ORIENTATION_UPSTREAM;
import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class LinxBreakendInterpreter
{
    @NotNull
    private final Map<String, StructuralVariant> structuralVariantsMap;
    @NotNull
    private final Map<Integer, LinxSvAnnotation> linxSvAnnotationsMap;
    @NotNull
    private final EnsemblDataCache ensemblDataCache;

    public LinxBreakendInterpreter(
            @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<LinxSvAnnotation> linxSvAnnotations,
            @NotNull final EnsemblDataCache ensemblDataCache)
    {
        this.structuralVariantsMap = structuralVariants.stream().collect(Collectors.toMap(StructuralVariant::id, s -> s));
        this.linxSvAnnotationsMap = linxSvAnnotations.stream().collect(Collectors.toMap(LinxSvAnnotation::svId, s -> s));
        this.ensemblDataCache = ensemblDataCache;
    }

    @NotNull
    public LinxBreakend interpret(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        LinxSvAnnotation svAnnotation = linxSvAnnotationsMap.get(linxBreakend.svId());
        if(svAnnotation == null)
        {
            LOGGER.warn("No linx sv annotation found for breakend {}", linxBreakend.id());
        }

        StructuralVariant structuralVariant = svAnnotation != null ? structuralVariantsMap.get(svAnnotation.vcfId()) : null;
        if(structuralVariant == null)
        {
            LOGGER.warn("No structural variant found for breakend {}", linxBreakend.id());
        }

        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .chromosome(chromosome(structuralVariant, linxBreakend.isStart()))
                .chromosomeBand(chromosomeBand(linxBreakend.gene()))
                .transcript(linxBreakend.transcriptId())
                .isCanonical(linxBreakend.canonical())
                .geneOrientation(linxBreakend.geneOrientation().equals(BREAKEND_ORIENTATION_UPSTREAM) ? LinxGeneOrientation.Upstream : LinxGeneOrientation.Downstream)
                .isCanonical(linxBreakend.canonical())
                .disruptive(linxBreakend.disruptive())
                .reported(linxBreakend.reportedDisruption())
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .type(svType(structuralVariant))
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .orientation(orientation(structuralVariant, linxBreakend.isStart()))
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                .junctionCopyNumber(junctionCopyNumber(svAnnotation))
                .build();
    }

    @NotNull
    private static String chromosome(@Nullable StructuralVariant structuralVariant, boolean isStart)
    {
        if(structuralVariant == null || structuralVariant.chromosome(isStart) == null)
        {
            return "";
        }
        else
        {
            String chromosome = structuralVariant.chromosome(isStart);
            return chromosome != null ? chromosome : "";
        }
    }

    @NotNull
    private String chromosomeBand(@NotNull String gene)
    {
        GeneData geneData = ensemblDataCache.getGeneDataByName(gene);
        return geneData != null ? geneData.KaryotypeBand : "";
    }

    @NotNull
    private static LinxBreakendType svType(@Nullable StructuralVariant structuralVariant)
    {
        if(structuralVariant == null)
        {
            return LinxBreakendType.BND;
        }
        else
        {
            return LinxBreakendType.valueOf(structuralVariant.type().name());
        }
    }

    private static byte orientation(@Nullable StructuralVariant structuralVariant, boolean isStart)
    {
        if(structuralVariant == null)
        {
            return 0;
        }
        else
        {
            Byte orientation = structuralVariant.orientation(isStart);
            return orientation != null ? orientation : 0;
        }
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