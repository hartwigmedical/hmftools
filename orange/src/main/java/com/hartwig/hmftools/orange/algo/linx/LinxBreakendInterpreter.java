package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.google.common.annotations.VisibleForTesting;
import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;

import org.jetbrains.annotations.NotNull;

public class LinxBreakendInterpreter
{

    @NotNull
    List<StructuralVariant> structuralVariants;
    @NotNull
    List<LinxSvAnnotation> linxSvAnnotations;
    @NotNull EnsemblDataCache ensemblDataCache;

    public LinxBreakendInterpreter(
            @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<LinxSvAnnotation> linxSvAnnotations,
            @NotNull final EnsemblDataCache ensemblDataCache
    )
    {
        this.structuralVariants = structuralVariants;
        this.linxSvAnnotations = linxSvAnnotations;
        this.ensemblDataCache = ensemblDataCache;
    }

    public LinxBreakend interpret(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        ImmutableLinxBreakend.Builder builder = breakendBuilder(linxBreakend);
        LinxSvAnnotation svAnnotation = linxSvAnnotations.stream()
                .filter(s -> s.svId() == linxBreakend.svId())
                .findFirst()
                .orElse(null);

        if(svAnnotation == null)
        {
            LOGGER.warn("No linx sv annotation found for breakend {}", linxBreakend.id());
            return builder.build();
        }

        StructuralVariant sv = structuralVariants.stream()
                .filter(s -> s.id().equals(svAnnotation.vcfId()))
                .findFirst()
                .orElse(null);

        if(sv == null)
        {
            LOGGER.warn("No structural variant found for breakend {}", linxBreakend.id());
            return builder.build();
        }

        String chrom = sv.chromosome(linxBreakend.isStart());
        if(chrom == null)
        {
            LOGGER.warn("No chromosome found for breakend {}", linxBreakend.id());
            return builder.build();
        }

        Byte orientation = sv.orientation(linxBreakend.isStart());
        if(orientation == null)
        {
            LOGGER.warn("No orientation found for breakend {}", linxBreakend.id());
            return builder.build();
        }

        double junctionCopyNumber = junctionCopyNumber(linxBreakend, sv, svAnnotation);
        GeneData geneData = ensemblDataCache.getGeneDataByName(linxBreakend.gene());
        String chrBand = geneData != null ? geneData.KaryotypeBand : "";

        return builder.chromosome(chrom)
                .chromosomeBand(chrBand)
                .type(LinxBreakendType.valueOf(sv.type().name()))
                .orientation(orientation)
                .junctionCopyNumber(junctionCopyNumber)
                .build();
    }

    @VisibleForTesting
    static double junctionCopyNumber(
            @NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend,
            @NotNull StructuralVariant sv,
            @NotNull LinxSvAnnotation svAnnotation
    )
    {
        if(sv.type() == StructuralVariantType.SGL && linxBreakend.regionType() == com.hartwig.hmftools.common.gene.TranscriptRegionType.IG)
        {
            return 0D;
        }
        else if(svAnnotation.junctionCopyNumberMin() == 0D)
        {
            return svAnnotation.junctionCopyNumberMax();
        }
        else
        {
            return 0.5D * (svAnnotation.junctionCopyNumberMin() + svAnnotation.junctionCopyNumberMax());
        }
    }

    @NotNull
    private static ImmutableLinxBreakend.Builder breakendBuilder(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .transcript(linxBreakend.transcriptId())
                .isCanonical(linxBreakend.canonical())
                .geneOrientation(linxBreakend.geneOrientation())
                .isCanonical(linxBreakend.canonical())
                .disruptive(linxBreakend.disruptive())
                .reported(linxBreakend.reportedDisruption())
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                // Placeholders values for fields to be derived from other sources, as fallback when derivation not possible
                .chromosome("")
                .chromosomeBand("")
                .type(LinxBreakendType.BND)
                .orientation(0)
                .junctionCopyNumber(0);
    }
}