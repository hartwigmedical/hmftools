package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxFusion;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.LinxSvAnnotation;

import org.jetbrains.annotations.NotNull;

public final class LinxConversion
{
    @NotNull
    public static LinxSvAnnotation convert(@NotNull com.hartwig.hmftools.common.linx.LinxSvAnnotation linxSvAnnotation)
    {
        return ImmutableLinxSvAnnotation.builder()
                .vcfId(linxSvAnnotation.vcfId())
                .svId(linxSvAnnotation.svId())
                .clusterId(linxSvAnnotation.clusterId())
                .clusterReason(linxSvAnnotation.clusterReason())
                .fragileSiteStart(linxSvAnnotation.fragileSiteStart())
                .fragileSiteEnd(linxSvAnnotation.fragileSiteEnd())
                .isFoldback(linxSvAnnotation.isFoldback())
                .lineTypeStart(linxSvAnnotation.lineTypeStart())
                .lineTypeEnd(linxSvAnnotation.lineTypeEnd())
                .junctionCopyNumberMin(linxSvAnnotation.junctionCopyNumberMin())
                .junctionCopyNumberMax(linxSvAnnotation.junctionCopyNumberMax())
                .geneStart(linxSvAnnotation.geneStart())
                .geneEnd(linxSvAnnotation.geneEnd())
                .localTopologyIdStart(linxSvAnnotation.localTopologyIdStart())
                .localTopologyIdEnd(linxSvAnnotation.localTopologyIdEnd())
                .localTopologyStart(linxSvAnnotation.localTopologyStart())
                .localTopologyEnd(linxSvAnnotation.localTopologyEnd())
                .localTICountStart(linxSvAnnotation.localTICountStart())
                .localTICountEnd(linxSvAnnotation.localTICountEnd())
                .build();
    }

    @NotNull
    public static LinxFusion convert(@NotNull com.hartwig.hmftools.common.linx.LinxFusion linxFusion)
    {
        return ImmutableLinxFusion.builder()
                .geneStart(linxFusion.geneStart())
                .geneContextStart(linxFusion.geneContextStart())
                .geneTranscriptStart(linxFusion.geneTranscriptStart())
                .geneEnd(linxFusion.geneEnd())
                .geneContextEnd(linxFusion.geneContextEnd())
                .geneTranscriptEnd(linxFusion.geneTranscriptEnd())
                .reported(linxFusion.reported())
                .reportedType(LinxFusionType.valueOf(linxFusion.reportedType()))
                .phased(FusionPhasedType.valueOf(linxFusion.phased().name()))
                .driverLikelihood(FusionLikelihoodType.valueOf(linxFusion.likelihood().name()))
                .fusedExonUp(linxFusion.fusedExonUp())
                .fusedExonDown(linxFusion.fusedExonDown())
                .chainLinks(linxFusion.chainLinks())
                .chainTerminated(linxFusion.chainTerminated())
                .domainsKept(linxFusion.domainsKept())
                .domainsLost(linxFusion.domainsLost())
                .junctionCopyNumber(linxFusion.junctionCopyNumber())
                .build();
    }

    @NotNull
    public static LinxBreakend convert(@NotNull com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend)
    {
        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .chromosome(linxBreakend.chromosome())
                .chromosomeBand(linxBreakend.chrBand())
                .transcript(linxBreakend.transcriptId())
                .isCanonical(linxBreakend.canonical())
                .geneOrientation(linxBreakend.geneOrientation())
                .isCanonical(linxBreakend.canonical())
                .orientation(linxBreakend.orientation())
                .disruptive(linxBreakend.disruptive())
                .reported(linxBreakend.reportedDisruption())
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .type(LinxBreakendType.valueOf(linxBreakend.type().name()))
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .orientation(linxBreakend.orientation())
                .strand(linxBreakend.strand())
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                .junctionCopyNumber(linxBreakend.junctionCopyNumber())
                .build();
    }

    @NotNull
    public static LinxHomozygousDisruption convert(@NotNull com.hartwig.hmftools.common.linx.HomozygousDisruption homozygousDisruption)
    {
        return ImmutableLinxHomozygousDisruption.builder()
                .gene(homozygousDisruption.gene())
                .chromosome(homozygousDisruption.chromosome())
                .chromosomeBand(homozygousDisruption.chromosomeBand())
                .transcript(homozygousDisruption.transcript())
                .isCanonical(homozygousDisruption.isCanonical())
                .build();
    }
}
