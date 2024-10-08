package com.hartwig.hmftools.orange.algo.linx;

import static com.hartwig.hmftools.orange.OrangeApplication.LOGGER;

import java.util.List;

import com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

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

    public LinxBreakend convert(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        LinxSvAnnotation svAnnotation = linxSvAnnotations.stream()
                .filter(s -> s.svId() == linxBreakend.svId())
                .findFirst()
                .orElse(null);

        if(svAnnotation == null)
        {
            LOGGER.warn("No linx sv annotation found for breakend {}", linxBreakend.id());
            return LinxConversion.convert(linxBreakend);
        }

        StructuralVariant sv = structuralVariants.stream()
                .filter(s -> s.id().equals(svAnnotation.vcfId()))
                .findFirst()
                .orElse(null);

        if(sv == null)
        {
            LOGGER.warn("No structural variant found for breakend {}", linxBreakend.id());
            return LinxConversion.convert(linxBreakend);
        }

        String chrom = sv.chromosome(linxBreakend.isStart());
        if(chrom == null)
        {
            LOGGER.warn("No chromosome found for breakend {}", linxBreakend.id());
            return LinxConversion.convert(linxBreakend);
        }

        Byte orientation = sv.orientation(linxBreakend.isStart());
        if(orientation == null)
        {
            LOGGER.warn("No orientation found for breakend {}", linxBreakend.id());
            return LinxConversion.convert(linxBreakend);
        }

        double junctionCopyNumber = junctionCopyNumber(linxBreakend, sv, svAnnotation);
        GeneData geneData = ensemblDataCache.getGeneDataByName(linxBreakend.gene());
        String chrBand = geneData != null ? geneData.KaryotypeBand : "";

        return ImmutableLinxBreakend.builder()
                .from(LinxConversion.convert(linxBreakend))
                .chromosome(chrom)
                .chromosomeBand(chrBand)
                .type(LinxBreakendType.valueOf(sv.type().name()))
                .orientation(orientation)
                .junctionCopyNumber(junctionCopyNumber)
                .build();
    }

    private static double junctionCopyNumber(
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
}