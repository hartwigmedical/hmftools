package com.hartwig.hmftools.orange.report;

import com.hartwig.hmftools.common.chord.ChordData;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.hla.LilacSummaryData;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.variant.AllelicDepth;
import com.hartwig.hmftools.datamodel.chord.ChordRecord;
import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.chord.ImmutableChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.linx.*;
import com.hartwig.hmftools.datamodel.orange.*;
import com.hartwig.hmftools.datamodel.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.datamodel.peach.ImmutablePeachRecord;
import com.hartwig.hmftools.datamodel.peach.PeachGenotype;
import com.hartwig.hmftools.datamodel.peach.PeachRecord;
import com.hartwig.hmftools.datamodel.purple.*;
import com.hartwig.hmftools.datamodel.sv.LinxBreakendType;
import com.hartwig.hmftools.datamodel.virus.*;
import com.hartwig.hmftools.orange.algo.OrangeReport;
import com.hartwig.hmftools.orange.algo.linx.LinxInterpretedData;
import com.hartwig.hmftools.orange.algo.purple.PurityPloidyFit;
import com.hartwig.hmftools.orange.algo.purple.PurpleInterpretedData;

import java.util.List;
import java.util.Objects;
import java.util.Optional;

public class OrangeReportToRecordConversion {

    private OrangeReportToRecordConversion() {
    }

    public static OrangeRecord convert(OrangeReport report) {
        RefGenomeVersion refGenomeVersion = report.refGenomeVersion();
        return ImmutableOrangeRecord.builder()
                .sampleId(report.sampleId())
                .experimentDate(report.experimentDate())
                .refGenomeVersion(OrangeRefGenomeVersion.valueOf(refGenomeVersion.name()))
                .purple(convert(report.purple()))
                .linx(convert(report.linx()))
                .lilac(convert(report.lilac()))
                .virusInterpreter(Optional.ofNullable(report.virusInterpreter()).map(OrangeReportToRecordConversion::convert).orElse(null))
                .chord(Optional.ofNullable(report.chord()).map(OrangeReportToRecordConversion::convert).orElse(null))
                .cuppa(Optional.ofNullable(report.cuppa()).map(OrangeReportToRecordConversion::convert).orElse(null))
                .peach(Optional.ofNullable(report.peach()).map(OrangeReportToRecordConversion::convert).orElse(null))
                .plots(convert(report.plots()))
                .build();
    }

    private static OrangePlots convert(com.hartwig.hmftools.orange.algo.OrangePlots plots) {
        return ImmutableOrangePlots.builder()
                .purpleFinalCircosPlot(plots.purpleFinalCircosPlot())
                .build();
    }

    private static PeachRecord convert(List<com.hartwig.hmftools.common.peach.PeachGenotype> peachGenotypes) {
        return ImmutablePeachRecord.builder()
                .entries(() -> peachGenotypes.stream().map(OrangeReportToRecordConversion::convert).iterator())
                .build();
    }

    private static PeachGenotype convert(com.hartwig.hmftools.common.peach.PeachGenotype peachGenotype) {
        return ImmutablePeachGenotype.builder()
                .gene(peachGenotype.gene())
                .haplotype(peachGenotype.haplotype())
                .function(peachGenotype.function())
                .linkedDrugs(peachGenotype.linkedDrugs())
                .urlPrescriptionInfo(peachGenotype.urlPrescriptionInfo())
                .panelVersion(peachGenotype.panelVersion())
                .repoVersion(peachGenotype.repoVersion())
                .build();
    }

    private static CuppaData convert(com.hartwig.hmftools.orange.algo.cuppa.CuppaData cuppaData) {
        return ImmutableCuppaData.builder()
                .predictions(() -> cuppaData.predictions().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .simpleDups32To200B(cuppaData.simpleDups32To200B())
                .maxComplexSize(cuppaData.maxComplexSize())
                .telomericSGLs(cuppaData.telomericSGLs())
                .LINECount(cuppaData.LINECount())
                .build();
    }

    private static CuppaPrediction convert(com.hartwig.hmftools.orange.algo.cuppa.CuppaPrediction cuppaPrediction) {
        return ImmutableCuppaPrediction.builder()
                .cancerType(cuppaPrediction.cancerType())
                .likelihood(cuppaPrediction.likelihood())
                .build();
    }

    private static ChordRecord convert(ChordData chordData) {
        return ImmutableChordRecord.builder()
                .hrdValue(chordData.hrdValue())
                .hrStatus(ChordStatus.valueOf(chordData.hrStatus().name()))
                .build();
    }

    private static VirusInterpreterData convert(com.hartwig.hmftools.common.virus.VirusInterpreterData interpreterData) {
        return ImmutableVirusInterpreterData.builder()
                .allViruses(() -> interpreterData.allViruses().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .reportableViruses(() -> interpreterData.reportableViruses().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .build();
    }

    private static AnnotatedVirus convert(com.hartwig.hmftools.common.virus.AnnotatedVirus annotatedVirus) {
        return ImmutableAnnotatedVirus.builder()
                .name(annotatedVirus.name())
                .qcStatus(VirusBreakendQCStatus.valueOf(annotatedVirus.qcStatus().name()))
                .integrations(annotatedVirus.integrations())
                .interpretation(annotatedVirus.interpretation())
                .percentageCovered(annotatedVirus.percentageCovered())
                .reported(annotatedVirus.reported())
                .build();
    }

    private static LilacRecord convert(LilacSummaryData lilacSummaryData) {
        return ImmutableLilacRecord.builder()
                .qc(lilacSummaryData.qc())
                .alleles(() -> lilacSummaryData.alleles().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .build();
    }

    private static LilacAllele convert(com.hartwig.hmftools.common.hla.LilacAllele allele) {
        return ImmutableLilacAllele.builder()
                .allele(allele.allele())
                .tumorCopyNumber(allele.tumorCopyNumber())
                .somaticMissense(allele.somaticMissense())
                .somaticNonsenseOrFrameshift(allele.somaticNonsenseOrFrameshift())
                .somaticSplice(allele.somaticSplice())
                .somaticSynonymous(allele.somaticSynonymous())
                .somaticInframeIndel(allele.somaticInframeIndel())
                .build();
    }

    private static LinxRecord convert(LinxInterpretedData linxInterpretedData) {
        return ImmutableLinxRecord.builder()
                .allSomaticStructuralVariants(() -> linxInterpretedData.allSomaticStructuralVariants().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .allSomaticFusions(() -> linxInterpretedData.allSomaticFusions().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .reportableSomaticFusions(() -> linxInterpretedData.reportableSomaticFusions().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .additionalSuspectSomaticFusions(() -> linxInterpretedData.additionalSuspectSomaticFusions().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .allSomaticBreakends(() -> linxInterpretedData.allSomaticBreakends().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .reportableSomaticBreakends(() -> linxInterpretedData.reportableSomaticBreakends().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .additionalSuspectSomaticBreakends(() -> linxInterpretedData.additionalSuspectSomaticBreakends().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .somaticHomozygousDisruptions(() -> linxInterpretedData.somaticHomozygousDisruptions().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .build();
    }

    private static LinxSvAnnotation convert(com.hartwig.hmftools.common.linx.LinxSvAnnotation linxSvAnnotation) {
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

    private static LinxFusion convert(com.hartwig.hmftools.common.linx.LinxFusion linxFusion) {
        return ImmutableLinxFusion.builder()
                .name(linxFusion.name())
                .reported(linxFusion.reported())
                .reportedType(linxFusion.reportedType())
                .phased(FusionPhasedType.valueOf(linxFusion.phased().name()))
                .likelihood(FusionLikelihoodType.valueOf(linxFusion.likelihood().name()))
                .fusedExonUp(linxFusion.fusedExonUp())
                .fusedExonDown(linxFusion.fusedExonDown())
                .geneStart(linxFusion.geneStart())
                .geneContextStart(linxFusion.geneContextStart())
                .geneTranscriptStart(linxFusion.geneTranscriptStart())
                .geneEnd(linxFusion.geneEnd())
                .geneContextEnd(linxFusion.geneContextEnd())
                .geneTranscriptEnd(linxFusion.geneTranscriptEnd())
                .junctionCopyNumber(linxFusion.junctionCopyNumber())
                .build();
    }

    private static LinxBreakend convert(com.hartwig.hmftools.common.linx.LinxBreakend linxBreakend) {
        return ImmutableLinxBreakend.builder()
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

    private static HomozygousDisruption convert(com.hartwig.hmftools.common.linx.HomozygousDisruption homozygousDisruption) {
        return ImmutableHomozygousDisruption.builder()
                .chromosome(homozygousDisruption.chromosome())
                .chromosomeBand(homozygousDisruption.chromosomeBand())
                .gene(homozygousDisruption.gene())
                .transcript(homozygousDisruption.transcript())
                .isCanonical(homozygousDisruption.isCanonical())
                .build();
    }

    private static PurpleRecord convert(PurpleInterpretedData purpleInterpretedData) {
        var somaticDriverIterator = purpleInterpretedData.somaticDrivers().stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        var germlineDriverIterator = Objects.requireNonNullElseGet(purpleInterpretedData.germlineDrivers(), List::<DriverCatalog>of).stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        var somaticVariantIterator = purpleInterpretedData.allSomaticVariants().stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        var germlineVariantIterator = Objects.requireNonNullElseGet(purpleInterpretedData.allGermlineVariants(), List::<com.hartwig.hmftools.orange.algo.purple.PurpleVariant>of).stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        var reportableGermlineVariantIterator = Objects.requireNonNullElseGet(purpleInterpretedData.reportableGermlineVariants(), List::<com.hartwig.hmftools.orange.algo.purple.PurpleVariant>of).stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        return ImmutablePurpleRecord.builder()
                .fit(convert(purpleInterpretedData.fit()))
                .characteristics(convert(purpleInterpretedData.characteristics()))
                .somaticDrivers(() -> somaticDriverIterator)
                .germlineDrivers(() -> germlineDriverIterator)
                .allSomaticVariants(() -> somaticVariantIterator)
                .reportableSomaticVariants(() -> purpleInterpretedData.reportableSomaticVariants().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .allGermlineVariants(() -> germlineVariantIterator)
                .reportableGermlineVariants(() -> reportableGermlineVariantIterator)
                .allSomaticCopyNumbers(() -> purpleInterpretedData.allSomaticCopyNumbers().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .allSomaticGeneCopyNumbers(() -> purpleInterpretedData.allSomaticGeneCopyNumbers().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .suspectGeneCopyNumbersWithLOH(() -> purpleInterpretedData.suspectGeneCopyNumbersWithLOH().stream().map(OrangeReportToRecordConversion::convert).iterator())
                .allSomaticGainsLosses(purpleInterpretedData.reportableSomaticGainsLosses())
                .reportableSomaticGainsLosses(purpleInterpretedData.reportableSomaticGainsLosses())
                .build();
    }

    private static PurpleFit convert(PurityPloidyFit fit) {
        var qcStatusIterator = fit.qc().status().stream()
                .map(status -> PurpleQCStatus.valueOf(status.name()))
                .iterator();
        return ImmutablePurpleFit.builder()
                .qcStatus(() -> qcStatusIterator)
                .hasSufficientQuality(fit.hasSufficientQuality())
                .containsTumorCells(fit.containsTumorCells())
                .purity(fit.purity())
                .ploidy(fit.ploidy())
                .build();
    }

    private static PurpleCharacteristics convert(com.hartwig.hmftools.orange.algo.purple.PurpleCharacteristics characteristics) {
        return ImmutablePurpleCharacteristics.builder()
                .microsatelliteIndelsPerMb(characteristics.microsatelliteIndelsPerMb())
                .microsatelliteStatus(PurpleMicrosatelliteStatus.valueOf(characteristics.microsatelliteStatus().name()))
                .tumorMutationalBurdenPerMb(characteristics.tumorMutationalBurdenPerMb())
                .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.valueOf(characteristics.tumorMutationalBurdenStatus().name()))
                .tumorMutationalLoad(characteristics.tumorMutationalLoad())
                .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.valueOf(characteristics.tumorMutationalLoadStatus().name()))
                .build();
    }

    private static PurpleDriver convert(DriverCatalog catalog) {
        return ImmutablePurpleDriver.builder()
                .gene(catalog.gene())
                .transcript(catalog.transcript())
                .driver(PurpleDriverType.valueOf(catalog.driver().name()))
                .transcript(catalog.transcript())
                .build();
    }

    private static PurpleVariant convert(com.hartwig.hmftools.orange.algo.purple.PurpleVariant variant) {
        var otherImpactsIterator = variant.otherImpacts().stream()
                .map(OrangeReportToRecordConversion::convert)
                .iterator();
        return ImmutablePurpleVariant.builder()
                .type(PurpleVariantType.valueOf(variant.type().name()))
                .gene(variant.gene())
                .chromosome(variant.chromosome())
                .position(variant.position())
                .ref(variant.ref())
                .alt(variant.alt())
                .canonicalImpact(convert(variant.canonicalImpact()))
                .otherImpacts(() -> otherImpactsIterator)
                .hotspot(Hotspot.valueOf(variant.hotspot().name()))
                .reported(variant.reported())
                .tumorDepth(convert(variant.tumorDepth()))
                .adjustedCopyNumber(variant.adjustedCopyNumber())
                .minorAlleleCopyNumber(variant.minorAlleleCopyNumber())
                .variantCopyNumber(variant.variantCopyNumber())
                .biallelic(variant.biallelic())
                .genotypeStatus(PurpleGenotypeStatus.valueOf(variant.genotypeStatus().name()))
                .subclonalLikelihood(variant.subclonalLikelihood())
                .localPhaseSets(variant.localPhaseSets())
                .build();
    }

    private static PurpleTranscriptImpact convert(com.hartwig.hmftools.orange.algo.purple.PurpleTranscriptImpact transcriptImpact) {
        var effectIterator = transcriptImpact.effects().stream().map(effect -> PurpleVariantEffect.valueOf(effect.name())).iterator();
        return ImmutablePurpleTranscriptImpact.builder()
                .transcript(transcriptImpact.transcript())
                .hgvsCodingImpact(transcriptImpact.hgvsCodingImpact())
                .hgvsProteinImpact(transcriptImpact.hgvsProteinImpact())
                .affectedCodon(transcriptImpact.affectedCodon())
                .affectedExon(transcriptImpact.affectedExon())
                .spliceRegion(transcriptImpact.spliceRegion())
                .effects(() -> effectIterator)
                .codingEffect(PurpleCodingEffect.valueOf(transcriptImpact.codingEffect().name()))
                .build();
    }

    private static PurpleAllelicDepth convert(AllelicDepth allelicDepth) {
        return ImmutablePurpleAllelicDepth.builder()
                .alleleReadCount(allelicDepth.alleleReadCount())
                .totalReadCount(allelicDepth.totalReadCount())
                .build();
    }

    private static PurpleCopyNumber convert(com.hartwig.hmftools.common.purple.PurpleCopyNumber copyNumber) {
        return ImmutablePurpleCopyNumber.builder()
                .chromosome(copyNumber.chromosome())
                .start(copyNumber.start())
                .end(copyNumber.end())
                .build();
    }

    private static PurpleGeneCopyNumber convert(GeneCopyNumber geneCopyNumber) {
        return ImmutablePurpleGeneCopyNumber.builder()
                .chromosome(geneCopyNumber.chromosome())
                .chromosomeBand(geneCopyNumber.chromosomeBand())
                .gene(geneCopyNumber.geneName())
                .minCopyNumber(geneCopyNumber.minCopyNumber())
                .minMinorAlleleCopyNumber(geneCopyNumber.minMinorAlleleCopyNumber())
                .build();
    }
}
