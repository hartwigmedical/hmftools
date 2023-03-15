package com.hartwig.hmftools.datamodel.linx;

import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestLinxFactory {

    private TestLinxFactory() {
    }

    @NotNull
    public static ImmutableLinxSvAnnotation.Builder structuralVariantBuilder() {
        return ImmutableLinxSvAnnotation.builder()
                .vcfId(Strings.EMPTY)
                .svId(0)
                .clusterId(0)
                .clusterReason(Strings.EMPTY)
                .fragileSiteStart(false)
                .fragileSiteEnd(false)
                .isFoldback(false)
                .lineTypeStart(Strings.EMPTY)
                .lineTypeEnd(Strings.EMPTY)
                .junctionCopyNumberMin(0)
                .junctionCopyNumberMax(0)
                .geneStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .localTopologyIdStart(0)
                .localTopologyIdEnd(0)
                .localTopologyStart(Strings.EMPTY)
                .localTopologyEnd(Strings.EMPTY)
                .localTICountStart(0)
                .localTICountEnd(0)
                ;
    }

    @NotNull
    public static ImmutableHomozygousDisruption.Builder homozygousDisruptionBuilder() {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false);
    }

    @NotNull
    public static ImmutableLinxBreakend.Builder breakendBuilder() {
        return ImmutableLinxBreakend.builder()
                .reportedDisruption(true)
                .disruptive(false)
                .id(0)
                .svId(0)
                .gene(Strings.EMPTY)
                .chromosome(Strings.EMPTY)
                .chrBand(Strings.EMPTY)
                .transcriptId(Strings.EMPTY)
                .canonical(false)
                .type(LinxBreakendType.BND)
                .junctionCopyNumber(0D)
                .undisruptedCopyNumber(0D)
                .nextSpliceExonRank(0)
                .exonUp(0)
                .exonDown(0)
                .geneOrientation(Strings.EMPTY)
                .orientation(0)
                .strand(0)
                .regionType(TranscriptRegionType.INTRONIC)
                .codingType(TranscriptCodingType.NON_CODING);
    }

    @NotNull
    public static ImmutableLinxFusion.Builder fusionBuilder() {
        return ImmutableLinxFusion.builder()
                .reported(true)
                .reportedType(LinxFusionType.NONE)
                .name(Strings.EMPTY)
                .geneStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .fusedExonUp(0)
                .geneEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .fusedExonDown(0)
                .likelihood(FusionLikelihoodType.LOW)
                .phased(FusionPhasedType.OUT_OF_FRAME)
                .junctionCopyNumber(0D)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY);
    }
}
