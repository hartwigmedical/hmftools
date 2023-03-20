package com.hartwig.hmftools.orange.conversion;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;

public class LinxConversion {

    private LinxConversion() {
    }

    public static LinxSvAnnotation convert(com.hartwig.hmftools.common.linx.LinxSvAnnotation linxSvAnnotation) {
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

    public static LinxFusion convert(com.hartwig.hmftools.common.linx.LinxFusion linxFusion) {
        return ImmutableLinxFusion.builder()
                .name(linxFusion.name())
                .reported(linxFusion.reported())
                .reportedType(LinxFusionType.valueOf(linxFusion.reportedType()))
                .phased(FusionPhasedType.valueOf(linxFusion.phased().name()))
                .likelihood(FusionLikelihoodType.valueOf(linxFusion.likelihood().name()))
                .fusedExonUp(linxFusion.fusedExonUp())
                .fusedExonDown(linxFusion.fusedExonDown())
                .chainLinks(linxFusion.chainLinks())
                .chainTerminated(linxFusion.chainTerminated())
                .domainsKept(linxFusion.domainsKept())
                .domainsLost(linxFusion.domainsLost())
                .geneStart(linxFusion.geneStart())
                .geneContextStart(linxFusion.geneContextStart())
                .geneTranscriptStart(linxFusion.geneTranscriptStart())
                .geneEnd(linxFusion.geneEnd())
                .geneContextEnd(linxFusion.geneContextEnd())
                .geneTranscriptEnd(linxFusion.geneTranscriptEnd())
                .junctionCopyNumber(linxFusion.junctionCopyNumber())
                .build();
    }

    public static LinxBreakend convert(com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend) {
        return ImmutableLinxBreakend.builder()
                .id(linxBreakend.id())
                .svId(linxBreakend.svId())
                .gene(linxBreakend.gene())
                .transcriptId(linxBreakend.transcriptId())
                .canonical(linxBreakend.canonical())
                .geneOrientation(linxBreakend.geneOrientation())
                .canonical(linxBreakend.canonical())
                .orientation(linxBreakend.orientation())
                .disruptive(linxBreakend.disruptive())
                .reportedDisruption(linxBreakend.reportedDisruption())
                .undisruptedCopyNumber(linxBreakend.undisruptedCopyNumber())
                .regionType(TranscriptRegionType.valueOf(linxBreakend.regionType().name()))
                .codingType(TranscriptCodingType.valueOf(linxBreakend.codingType().name()))
                .nextSpliceExonRank(linxBreakend.nextSpliceExonRank())
                .type(LinxBreakendType.valueOf(linxBreakend.type().name()))
                .chromosome(linxBreakend.chromosome())
                .orientation(linxBreakend.orientation())
                .strand(linxBreakend.strand())
                .chrBand(linxBreakend.chrBand())
                .exonUp(linxBreakend.exonUp())
                .exonDown(linxBreakend.exonDown())
                .junctionCopyNumber(linxBreakend.junctionCopyNumber())
                .build();
    }

    public static HomozygousDisruption convert(com.hartwig.hmftools.common.linx.HomozygousDisruption homozygousDisruption) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(homozygousDisruption.chromosome())
                .chromosomeBand(homozygousDisruption.chromosomeBand())
                .gene(homozygousDisruption.gene())
                .transcript(homozygousDisruption.transcript())
                .isCanonical(homozygousDisruption.isCanonical())
                .build();
    }
}
