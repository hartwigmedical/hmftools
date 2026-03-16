package com.hartwig.hmftools.orange.algo.linx;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.common.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.common.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.common.linx.LinxTestFactory;

public final class TestLinxFactory
{
    public static ImmutableLinxFusion.Builder fusionBuilder()
    {
        return LinxTestFactory.fusionBuilder();

        /*
        return ImmutableLinxFusion.builder()
                .fivePrimeBreakendId(0)
                .threePrimeBreakendId(0)
                .name("GENE01_GENE02")
                .reported(true)
                .reportedType(KnownFusionType.EXON_DEL_DUP.toString())
                .reportableReasons(Collections.emptyList())
                .phased(FusionPhasedType.INFRAME)
                .likelihood(FusionLikelihoodType.HIGH)
                .fivePrimeVcfId("")
                .threePrimeVcfId("")
                .fivePrimeCoords("")
                .threePrimeCoords("")
                .chainLength(0)
                .chainLinks(0)
                .chainTerminated(false)
                .domainsKept("")
                .domainsLost("")
                .fusedExonUp(1)
                .fusedExonDown(2)
                .skippedExonsUp(0)
                .skippedExonsDown(0);
        */
    }

    public static ImmutableLinxBreakend.Builder breakendBuilder()
    {
        return LinxTestFactory.breakendBuilder();
        /*
        return ImmutableLinxBreakend.builder()
                .id(0)
                .svId(0)
                .vcfId("")
                .coords("")
                .gene("GENE")
                .isStart(true)
                .transcriptId(Strings.EMPTY)
                .canonical(true)
                .geneOrientation("")
                .disruptive(false)
                .reportedStatus(ReportedStatus.NONE)
                .undisruptedCopyNumber(0)
                .regionType(TranscriptRegionType.UNKNOWN)
                .codingType(TranscriptCodingType.UNKNOWN)
                .biotype("")
                .nextSpliceExonRank(0)
                .nextSpliceExonPhase(0)
                .nextSpliceDistance(0)
                .exonUp(0)
                .exonDown(0)
                .totalExonCount(10);
        */
    }

    public static ImmutableLinxSvAnnotation.Builder svAnnotationBuilder()
    {
        return LinxTestFactory.svAnnotationBuilder();

        /*
        return ImmutableLinxSvAnnotation.builder()
                .svId(1)
                .vcfIdStart("")
                .vcfIdEnd("")
                .coordsStart("")
                .coordsEnd("")
                .type(StructuralVariantType.BND)
                .clusterId(0)
                .clusterReason("")
                .fragileSiteStart(false)
                .fragileSiteEnd(false)
                .isFoldback(false)
                .lineTypeStart("")
                .lineTypeEnd("")
                .junctionCopyNumberMin(1.0D)
                .junctionCopyNumberMax(2.0D)
                .geneStart("")
                .geneEnd("")
                .localTopologyIdStart(0)
                .localTopologyIdEnd(0)
                .localTopologyStart("")
                .localTopologyEnd("")
                .localTICountStart(0)
                .localTICountEnd(0);
        */
    }

    public static ImmutableLinxData.Builder linxDataBuilder()
    {
        return ImmutableLinxData.builder()
                .somaticDrivers(Lists.newArrayList())
                .somaticBreakends(Lists.newArrayList())
                .somaticSvAnnotations(Lists.newArrayList())
                .somaticHomozygousDisruptions(Lists.newArrayList())
                .somaticDriverData(Lists.newArrayList())
                .somaticBreakendClusterIds(Maps.newHashMap())
                .fusions(Lists.newArrayList())
                .fusionClusterIds(Maps.newHashMap())
                .germlineDisruptions(Lists.newArrayList())
                .germlineBreakends(Lists.newArrayList())
                .germlineDrivers(Lists.newArrayList())
                .germlineSvAnnotations(Lists.newArrayList())
                .reportableEventPlots(Lists.newArrayList());
    }
}
