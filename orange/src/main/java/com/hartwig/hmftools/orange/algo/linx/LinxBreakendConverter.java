package com.hartwig.hmftools.orange.algo.linx;

import java.util.List;
import java.util.Objects;

import com.hartwig.hmftools.common.genome.chromosome.CytoBands;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantType;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.orange.conversion.LinxConversion;

import org.jetbrains.annotations.NotNull;

public class LinxBreakendConverter
{
    @NotNull
    List<StructuralVariant> structuralVariants;
    @NotNull
    List<LinxSvAnnotation> linxSvAnnotations;
    @NotNull
    CytoBands cytoBands;

    public LinxBreakendConverter(
            @NotNull final List<StructuralVariant> structuralVariants,
            @NotNull final List<LinxSvAnnotation> linxSvAnnotations,
            @NotNull RefGenomeVersion refGenomeVersion
    )
    {
        this.structuralVariants = structuralVariants;
        this.linxSvAnnotations = linxSvAnnotations;
        this.cytoBands = new CytoBands(refGenomeVersion);
    }

    public LinxBreakend convert(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        LinxSvAnnotation svAnnotation = linxSvAnnotations.stream()
                .filter(s -> s.svId() == linxBreakend.svId())
                .findFirst()
                .orElseThrow(() -> new RuntimeException("No sv annotation found for breakend " + linxBreakend.id()));

        if(svAnnotation.vcfId().startsWith("purple"))
        {
            System.out.println("  TODO: Skipping annotation vcfId that fails to map: " + svAnnotation.vcfId());
            return LinxConversion.convert(linxBreakend);
        }

        StructuralVariant sv = structuralVariants.stream()
                .filter(s -> s.id().equals(svAnnotation.vcfId()))
                .findFirst()
                .orElseThrow(() -> new RuntimeException("No structural variant found for breakend " + linxBreakend.id()));

        String chrom = sv.chromosome(linxBreakend.isStart());
        int position = Objects.requireNonNull(sv.position(linxBreakend.isStart()), "Position is null for breakend " + linxBreakend.id());
        byte orientation =
                Objects.requireNonNull(sv.orientation(linxBreakend.isStart()), "Orientation is null for breakend " + linxBreakend.id());

        StructuralVariantType type = sv.type();

        double junctionCopyNumber = 0.5D * (svAnnotation.junctionCopyNumberMin() + svAnnotation.junctionCopyNumberMax());
        String chrBand = cytoBands.getCytoBandName("chr" + chrom, position);

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
                .chromosome(chrom)
                .chromosomeBand(chrBand)
                .type(LinxBreakendType.valueOf(type.name()))
                .orientation(orientation)
                .junctionCopyNumber(junctionCopyNumber)
                .build();
    }
}