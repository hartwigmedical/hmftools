package com.hartwig.hmftools.orange.algo.linx;

import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.fusion.KnownFusionCache;
import com.hartwig.hmftools.common.linx.LinxData;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Objects;
import java.util.function.Supplier;

public class LinxInterpreter {

    private static final Logger LOGGER = LogManager.getLogger(LinxInterpreter.class);

    @NotNull
    private final List<DriverGene> driverGenes;
    @NotNull
    private final KnownFusionCache knownFusionCache;

    public LinxInterpreter(@NotNull final List<DriverGene> driverGenes, @NotNull final KnownFusionCache knownFusionCache) {
        this.driverGenes = driverGenes;
        this.knownFusionCache = knownFusionCache;
    }

    @NotNull
    public LinxRecord interpret(@NotNull LinxData linx) {
        LOGGER.info("Analysing linx data");
        List<com.hartwig.hmftools.common.linx.LinxFusion> additionalSuspectSomaticFusions =
                DNAFusionSelector.selectInterestingUnreportedFusions(linx.allSomaticFusions(), driverGenes);
        LOGGER.info(" Found an additional {} suspect somatic fusions that are potentially interesting",
                additionalSuspectSomaticFusions.size());

        List<com.hartwig.hmftools.common.linx.LinxBreakend> additionalSuspectSomaticBreakends =
                BreakendSelector.selectInterestingUnreportedBreakends(linx.allSomaticBreakends(),
                        linx.reportableSomaticFusions(),
                        knownFusionCache);
        LOGGER.info(" Found an additional {} suspect somatic breakends that are potentially interesting",
                additionalSuspectSomaticBreakends.size());

        return ImmutableLinxRecord.builder()
                .allSomaticStructuralVariants(() -> linx.allSomaticStructuralVariants().stream().map(LinxInterpreter::convert).iterator())
                .allGermlineStructuralVariants(() -> argumentOrEmpty(linx.allGermlineStructuralVariants()).stream().map(LinxInterpreter::convert).iterator())
                .allSomaticFusions(() -> linx.allSomaticFusions().stream().map(LinxInterpreter::convert).iterator())
                .reportableSomaticFusions(() -> linx.reportableSomaticFusions().stream().map(LinxInterpreter::convert).iterator())
                .additionalSuspectSomaticFusions(() -> additionalSuspectSomaticFusions.stream().map(LinxInterpreter::convert).iterator())
                .allSomaticBreakends(() -> linx.allSomaticBreakends().stream().map(LinxInterpreter::convert).iterator())
                .allGermlineBreakends(() -> argumentOrEmpty(linx.allGermlineBreakends()).stream().map(LinxInterpreter::convert).iterator())
                .reportableSomaticBreakends(() -> linx.reportableSomaticBreakends().stream().map(LinxInterpreter::convert).iterator())
                .reportableGermlineBreakends(() -> argumentOrEmpty(linx.reportableGermlineBreakends()).stream().map(LinxInterpreter::convert).iterator())
                .additionalSuspectSomaticBreakends(() -> additionalSuspectSomaticBreakends.stream().map(LinxInterpreter::convert).iterator())
                .somaticHomozygousDisruptions(() -> linx.somaticHomozygousDisruptions().stream().map(LinxInterpreter::convert).iterator())
                .build();
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
                .reportedType(linxFusion.reportedType())
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

    private static <T> List<T> argumentOrEmpty(List<T> obj) {
        return Objects.requireNonNullElseGet(obj, List::of);
    }
}
