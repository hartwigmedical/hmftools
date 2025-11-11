package com.hartwig.hmftools.common.linx;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.gene.TranscriptCodingType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;
import com.hartwig.hmftools.common.sv.StructuralVariantType;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LinxTestFactory
{
    @NotNull
    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder()
    {
        return ImmutableLinxSvAnnotation.builder()
                .vcfIdStart(Strings.EMPTY)
                .vcfIdEnd(Strings.EMPTY)
                .svId(0)
                .coordsStart(Strings.EMPTY)
                .coordsEnd(Strings.EMPTY)
                .clusterId(0)
                .clusterReason(Strings.EMPTY)
                .fragileSiteStart(false)
                .fragileSiteEnd(false)
                .isFoldback(false)
                .lineTypeStart(Strings.EMPTY)
                .lineTypeEnd(Strings.EMPTY)
                .junctionCopyNumberMin(0D)
                .junctionCopyNumberMax(0D)
                .geneStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .localTopologyIdStart(0)
                .localTopologyIdEnd(0)
                .localTopologyStart(Strings.EMPTY)
                .localTopologyEnd(Strings.EMPTY)
                .localTICountStart(0)
                .localTICountEnd(0);
    }

    @NotNull
    public static LinxFusion createMinimalTestFusion()
    {
        return fusionBuilder().build();
    }

    @NotNull
    public static ImmutableLinxFusion.Builder fusionBuilder()
    {
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(0)
                .threePrimeBreakendId(0)
                .fivePrimeVcfId(Strings.EMPTY)
                .threePrimeVcfId(Strings.EMPTY)
                .fivePrimeCoords(Strings.EMPTY)
                .threePrimeCoords(Strings.EMPTY)
                .name(Strings.EMPTY)
                .reported(false)
                .reportedType(KnownFusionType.NONE.toString())
                .reportableReasons("OK")
                .phased(FusionPhasedType.OUT_OF_FRAME)
                .likelihood(FusionLikelihoodType.NA)
                .chainLength(0)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept(Strings.EMPTY)
                .domainsLost(Strings.EMPTY)
                .skippedExonsUp(0)
                .skippedExonsDown(0)
                .fusedExonUp(0)
                .fusedExonDown(0)
                .geneStart(Strings.EMPTY)
                .geneContextStart(Strings.EMPTY)
                .geneTranscriptStart(Strings.EMPTY)
                .geneEnd(Strings.EMPTY)
                .geneContextEnd(Strings.EMPTY)
                .geneTranscriptEnd(Strings.EMPTY)
                .junctionCopyNumber(0D);
    }

    @NotNull
    public static ImmutableLinxBreakend.Builder breakendBuilder()
    {
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .vcfId(Strings.EMPTY)
                .coords(Strings.EMPTY)
                .isStart(true)
                .gene(Strings.EMPTY)
                .transcriptId(Strings.EMPTY)
                .canonical(true)
                .geneOrientation(Strings.EMPTY)
                .disruptive(false)
                .reportedDisruption(false)
                .undisruptedCopyNumber(0D)
                .regionType(TranscriptRegionType.UNKNOWN)
                .codingType(TranscriptCodingType.UNKNOWN)
                .biotype(Strings.EMPTY)
                .exonicBasePhase(0)
                .nextSpliceExonRank(0)
                .nextSpliceExonPhase(0)
                .nextSpliceDistance(0)
                .totalExonCount(0)
                .exonUp(0)
                .exonDown(0);
    }

    @NotNull
    public static ImmutableHomozygousDisruption.Builder homozygousDisruptionBuilder()
    {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(Strings.EMPTY)
                .chromosomeBand(Strings.EMPTY)
                .gene(Strings.EMPTY)
                .transcript(Strings.EMPTY)
                .isCanonical(false);
    }
}
