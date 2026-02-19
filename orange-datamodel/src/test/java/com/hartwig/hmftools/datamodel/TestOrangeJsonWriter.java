package com.hartwig.hmftools.datamodel;

import java.io.IOException;
import java.time.LocalDate;
import java.util.List;

import com.hartwig.hmftools.datamodel.chord.ChordStatus;
import com.hartwig.hmftools.datamodel.chord.ImmutableChordRecord;
import com.hartwig.hmftools.datamodel.cuppa.CuppaData;
import com.hartwig.hmftools.datamodel.cuppa.CuppaMode;
import com.hartwig.hmftools.datamodel.cuppa.CuppaPrediction;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaData;
import com.hartwig.hmftools.datamodel.cuppa.ImmutableCuppaPrediction;
import com.hartwig.hmftools.datamodel.driver.DriverInterpretation;
import com.hartwig.hmftools.datamodel.driver.ReportedStatus;
import com.hartwig.hmftools.datamodel.flagstat.ImmutableFlagstat;
import com.hartwig.hmftools.datamodel.gene.TranscriptCodingType;
import com.hartwig.hmftools.datamodel.gene.TranscriptRegionType;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.datamodel.hla.ImmutableLilacRecord;
import com.hartwig.hmftools.datamodel.hla.LilacAllele;
import com.hartwig.hmftools.datamodel.hla.LilacRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.immuno.ImmutableImmuneEscapeRecord;
import com.hartwig.hmftools.datamodel.linx.FusionLikelihoodType;
import com.hartwig.hmftools.datamodel.linx.FusionPhasedType;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxBreakend;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxFusion;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxHomozygousDisruption;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxRecord;
import com.hartwig.hmftools.datamodel.linx.ImmutableLinxSvAnnotation;
import com.hartwig.hmftools.datamodel.linx.LinxBreakend;
import com.hartwig.hmftools.datamodel.linx.LinxBreakendType;
import com.hartwig.hmftools.datamodel.linx.LinxFusionType;
import com.hartwig.hmftools.datamodel.linx.LinxGeneOrientation;
import com.hartwig.hmftools.datamodel.linx.LinxRecord;
import com.hartwig.hmftools.datamodel.linx.LinxUnreportableReason;
import com.hartwig.hmftools.datamodel.metrics.ImmutableWGSMetrics;
import com.hartwig.hmftools.datamodel.orange.ExperimentType;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangePlots;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeRecord;
import com.hartwig.hmftools.datamodel.orange.ImmutableOrangeSample;
import com.hartwig.hmftools.datamodel.orange.OrangeRecord;
import com.hartwig.hmftools.datamodel.orange.OrangeRefGenomeVersion;
import com.hartwig.hmftools.datamodel.orange.OrangeSample;
import com.hartwig.hmftools.datamodel.peach.ImmutablePeachGenotype;
import com.hartwig.hmftools.datamodel.purple.CopyNumberInterpretation;
import com.hartwig.hmftools.datamodel.purple.HotspotType;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleAllelicDepth;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCharacteristics;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleDriver;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleFit;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGainDeletion;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleGeneCopyNumber;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleQC;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleRecord;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleTranscriptImpact;
import com.hartwig.hmftools.datamodel.purple.ImmutablePurpleVariant;
import com.hartwig.hmftools.datamodel.purple.ImmutableTumorStats;
import com.hartwig.hmftools.datamodel.purple.PurpleCodingEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleDriver;
import com.hartwig.hmftools.datamodel.purple.PurpleDriverType;
import com.hartwig.hmftools.datamodel.purple.PurpleFittedPurityMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleGenotypeStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleLikelihoodMethod;
import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleQCStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleRecord;
import com.hartwig.hmftools.datamodel.purple.PurpleSomaticLikelihood;
import com.hartwig.hmftools.datamodel.purple.PurpleTumorMutationalStatus;
import com.hartwig.hmftools.datamodel.purple.PurpleVariant;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantEffect;
import com.hartwig.hmftools.datamodel.purple.PurpleVariantType;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.ImmutableVirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusBreakendQCStatus;
import com.hartwig.hmftools.datamodel.virus.VirusInterpretation;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterData;
import com.hartwig.hmftools.datamodel.virus.VirusInterpreterEntry;
import com.hartwig.hmftools.datamodel.virus.VirusLikelihoodType;

import org.jetbrains.annotations.NotNull;

// use this class to write the resources/minimally.populated.orange.json
public class TestOrangeJsonWriter
{
    public static void main(String[] args) throws IOException
    {
        // change the path to suit your set up
        String path = args[0];
        OrangeRecord record = createOrangeRecord();
        OrangeJson.getInstance().write(record, path);
    }

    @NotNull
    public static OrangeRecord createOrangeRecord()
    {
        return ImmutableOrangeRecord.builder()
                .sampleId("TEST")
                .samplingDate(LocalDate.of(2022, 1, 20))
                .experimentType(ExperimentType.WHOLE_GENOME)
                .refGenomeVersion(OrangeRefGenomeVersion.V37)
                .tumorSample(createOrangeSample())
                .purple(createPurpleRecord())
                .linx(createLinxRecord())
                .lilac(createLilacRecord())
                .immuneEscape(createImmuneEscapeRecord())
                .virusInterpreter(createVirusInterpreterData())
                .chord(ImmutableChordRecord.builder()
                        .brca1Value(0)
                        .brca2Value(0)
                        .hrdValue(0)
                        .hrdType("")
                        .hrStatus(ChordStatus.HR_PROFICIENT)
                        .build())
                .addPeach(ImmutablePeachGenotype.builder()
                        .gene("DPYD")
                        .allele("*1")
                        .alleleCount(2)
                        .haplotype("*1_HOM")
                        .function("Normal Function")
                        .linkedDrugs("5-Fluorouracil")
                        .urlPrescriptionInfo("https://www.pharmgkb.org/guidelineAnnotation/PA166104939")
                        .panelVersion("peach_prod_v1.3")
                        .repoVersion("1.7")
                        .build())
                .cuppa(createCuppaData())
                .plots(ImmutableOrangePlots.builder()
                        .purpleFinalCircosPlot("plot/empty.circos.png")
                        .sageTumorBQRPlot("")
                        .purpleInputPlot("")
                        .purpleClonalityPlot("")
                        .purpleCopyNumberPlot("")
                        .purpleVariantCopyNumberPlot("")
                        .purplePurityRangePlot("")
                        .purpleKataegisPlot("")
                        .build())
                .build();
    }

    @NotNull
    private static OrangeSample createOrangeSample() {
        return ImmutableOrangeSample.builder()
                .metrics(ImmutableWGSMetrics.builder()
                        .meanCoverage(0)
                        .sdCoverage(0)
                        .medianCoverage(0)
                        .madCoverage(0)
                        .pctExcAdapter(0.0)
                        .pctExcMapQ(0)
                        .pctExcDupe(0)
                        .pctExcUnpaired(0)
                        .pctExcBaseQ(0)
                        .pctExcOverlap(0)
                        .pctExcCapped(0)
                        .pctExcTotal(0)
                        .build())
                .flagstat(ImmutableFlagstat.builder()
                        .uniqueReadCount(0)
                        .secondaryCount(0)
                        .supplementaryCount(0)
                        .mappedProportion(0)
                        .build())
                .build();
    }

    @NotNull
    private static PurpleRecord createPurpleRecord()
    {
        PurpleDriver mutationDriver = ImmutablePurpleDriver.builder()
                        .gene("SF3B1")
                        .transcript("ENST00000335508")
                        .type(PurpleDriverType.MUTATION)
                        .driverLikelihood(0.2)
                        .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.HIGH)
                        .isCanonical(false)
                        .build();

        PurpleDriver deletionDriver = ImmutablePurpleDriver.builder()
                        .gene("SMAD4")
                        .transcript("ENST00000342988")
                        .type(PurpleDriverType.DEL)
                        .driverLikelihood(1.0)
                        .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.HIGH)
                        .isCanonical(false)
                        .build();

        PurpleVariant somaticVariant = ImmutablePurpleVariant.builder()
                .gene("SF3B1")
                .type(PurpleVariantType.SNP)
                .chromosome("2")
                .position(198266779)
                .ref("G")
                .alt("A")
                .worstCodingEffect(PurpleCodingEffect.MISSENSE)
                .variantCopyNumber(2.03)
                .minorAlleleCopyNumber(0.4)
                .adjustedCopyNumber(3.02)
                .adjustedVAF(1.2)
                .hotspot(HotspotType.NON_HOTSPOT)
                .tumorDepth(ImmutablePurpleAllelicDepth.builder().totalReadCount(20).alleleReadCount(10).build())
                .subclonalLikelihood(0.0)
                .somaticLikelihood(PurpleSomaticLikelihood.MEDIUM)
                .biallelic(false)
                .biallelicProbability(0.1)
                .genotypeStatus(PurpleGenotypeStatus.UNKNOWN)
                .repeatCount(1)
                .canonicalImpact(ImmutablePurpleTranscriptImpact.builder()
                        .transcript("ENST00000335508")
                        .hgvsCodingImpact("c.2153C>T")
                        .hgvsProteinImpact("p.Pro718Leu")
                        .affectedCodon(2153)
                        .affectedExon(12)
                        .inSpliceRegion(false)
                        .addEffects(PurpleVariantEffect.MISSENSE)
                        .codingEffect(PurpleCodingEffect.MISSENSE)
                        .reported(true)
                        .build())
                .build();

        PurpleVariant germlineVariant = ImmutablePurpleVariant.builder()
                .gene("BRCA1")
                .type(PurpleVariantType.SNP)
                .chromosome("17")
                .position(41209068)
                .ref("C")
                .alt("T")
                .worstCodingEffect(PurpleCodingEffect.SPLICE)
                .variantCopyNumber(1.0)
                .minorAlleleCopyNumber(0.8)
                .adjustedCopyNumber(2.0)
                .adjustedVAF(1.2)
                .hotspot(HotspotType.HOTSPOT)
                .tumorDepth(ImmutablePurpleAllelicDepth.builder().totalReadCount(30).alleleReadCount(20).build())
                .subclonalLikelihood(0.2)
                .somaticLikelihood(PurpleSomaticLikelihood.UNKNOWN)
                .biallelic(false)
                .biallelicProbability(0.0)
                .genotypeStatus(PurpleGenotypeStatus.HET)
                .repeatCount(1)
                .addLocalPhaseSets(1)
                .addLocalPhaseSets(2)
                .canonicalImpact(ImmutablePurpleTranscriptImpact.builder()
                        .transcript("ENST00000471181")
                        .hgvsCodingImpact("c.5340+1G>A")
                        .hgvsProteinImpact("p.?")
                        .inSpliceRegion(true)
                        .addEffects(PurpleVariantEffect.SPLICE_DONOR)
                        .addEffects(PurpleVariantEffect.INTRONIC)
                        .codingEffect(PurpleCodingEffect.SPLICE)
                        .reported(true)
                        .build())
                .build();


        return ImmutablePurpleRecord.builder()
                .fit(ImmutablePurpleFit.builder()
                        .qc(ImmutablePurpleQC.builder()
                                .addStatus(PurpleQCStatus.PASS)
                                .amberMeanDepth(111)
                                .contamination(0.0)
                                .totalCopyNumberSegments(100)
                                .unsupportedCopyNumberSegments(0)
                                .deletedGenes(4)
                                .build())
                        .fittedPurityMethod(PurpleFittedPurityMethod.NORMAL)
                        .purity(0.12)
                        .minPurity(0.97)
                        .maxPurity(1.0)
                        .ploidy(3.1)
                        .minPloidy(3.1)
                        .maxPloidy(3.15)
                        .build())
                .tumorStats(ImmutableTumorStats.builder()
                        .maxDiploidProportion(0.0211)
                        .hotspotMutationCount(3)
                        .hotspotStructuralVariantCount(0)
                        .smallVariantAlleleReadCount(2_273_196)
                        .structuralVariantTumorFragmentCount(4_908)
                        .bafCount(675_344)
                        .build())
                .characteristics(ImmutablePurpleCharacteristics.builder()
                        .wholeGenomeDuplication(true)
                        .microsatelliteIndelsPerMb(0.1)
                        .microsatelliteStatus(PurpleMicrosatelliteStatus.MSS)
                        .tumorMutationalBurdenPerMb(13.71)
                        .tumorMutationalBurdenStatus(PurpleTumorMutationalStatus.HIGH)
                        .tumorMutationalLoad(185)
                        .tumorMutationalLoadStatus(PurpleTumorMutationalStatus.HIGH)
                        .svTumorMutationalBurden(75)
                        .build())
                .somaticDrivers(List.of(mutationDriver, deletionDriver))
                .germlineDrivers(List.of(ImmutablePurpleDriver.builder()
                        .gene("BRCA1")
                        .transcript("ENST00000471181")
                        .type(PurpleDriverType.GERMLINE_MUTATION)
                        .driverLikelihood(0.8)
                        .likelihoodMethod(PurpleLikelihoodMethod.HOTSPOT)
                        .reportedStatus(ReportedStatus.REPORTED)
                        .driverInterpretation(DriverInterpretation.LOW)
                        .isCanonical(false)
                        .build()))
                .somaticVariants(List.of(somaticVariant))
                .germlineVariants(List.of(germlineVariant))
                .somaticCopyNumbers(List.of(ImmutablePurpleCopyNumber.builder()
                        .chromosome("1")
                        .start(10)
                        .end(20)
                        .averageTumorCopyNumber(4.1)
                        .build())
                )
                .somaticGeneCopyNumbers(List.of(ImmutablePurpleGeneCopyNumber.builder()
                        .gene("gene")
                        .chromosome("12")
                        .chromosomeBand("p13")
                        .transcript("trans")
                        .isCanonical(true)
                        .minCopyNumber(1.2)
                        .maxCopyNumber(1.2)
                        .minMinorAlleleCopyNumber(0.4)
                        .build())
                )
                .somaticGainsDels(List.of(ImmutablePurpleGainDeletion.builder()
                                .driver(deletionDriver)
                        .interpretation(CopyNumberInterpretation.FULL_DEL)
                        .chromosome("5")
                        .chromosomeBand("q2.2")
                        .minCopies(0.1)
                        .maxCopies(1.2)
                        .minMinorAlleleCopies(0.1)
                        .build()))
                .build();
    }

    @NotNull
    private static LinxRecord createLinxRecord()
    {
        LinxBreakend somaticBreakend = ImmutableLinxBreakend.builder()
                .reportedStatus(ReportedStatus.NOT_REPORTED)
                .disruptive(false)
                .id(0)
                .svId(1)
                .gene("NF1")
                .chromosome("1")
                .chromosomeBand("p12")
                .transcript("trans")
                .isCanonical(true)
                .type(LinxBreakendType.DUP)   // or: PurpleSVType.DUP — adjust if needed
                .junctionCopyNumber(1.1)
                .undisruptedCopyNumber(1.0)
                .nextSpliceExonRank(-1)
                .exonUp(1)
                .exonDown(2)
                .geneOrientation(LinxGeneOrientation.UPSTREAM)
                .orientation(-1)
                .regionType(TranscriptRegionType.EXONIC)
                .codingType(TranscriptCodingType.UTR_3P)
                .build();

        return ImmutableLinxRecord.builder()
                .addSomaticStructuralVariants(ImmutableLinxSvAnnotation.builder()
                        .vcfId("id")
                        .svId(1)
                        .clusterId(2)
                        .clusterReason("")
                        .fragileSiteStart(false)
                        .fragileSiteEnd(false)
                        .isFoldback(false)
                        .lineTypeStart("NONE")
                        .lineTypeEnd("NONE")
                        .junctionCopyNumberMin(2.0)
                        .junctionCopyNumberMax(3.0)
                        .geneStart("PTENR")
                        .geneEnd("PTEN")
                        .localTopologyIdStart(0)
                        .localTopologyIdEnd(1)
                        .localTopologyStart("ISOLATED_S")
                        .localTopologyEnd("ISOLATED_BE")
                        .localTICountStart(3)
                        .localTICountEnd(4)
                        .build())
                .addFusions(ImmutableLinxFusion.builder()
                        .reported(true)
                        .driverLikelihood(FusionLikelihoodType.HIGH)
                        .reportedType(LinxFusionType.KNOWN_PAIR)
                        .addUnreportedReasons(LinxUnreportableReason.NONE)
                        .fusedExonUp(1)
                        .fusedExonDown(2)
                        .geneStart("TMPRSS2")
                        .geneTranscriptStart("ENST00000332149")
                        .geneContextStart("Exon 1")
                        .geneEnd("ETV4")
                        .geneTranscriptEnd("ENST00000319349")
                        .geneContextEnd("Exon 2")
                        .phased(FusionPhasedType.INFRAME)
                        .junctionCopyNumber(1.1)
                        .chainLinks(0)
                        .chainTerminated(false)
                        .domainsKept("")
                        .domainsLost("")
                        .build())
                .addSomaticBreakends(somaticBreakend)
                .addSomaticHomozygousDisruptions(ImmutableLinxHomozygousDisruption.builder()
                        .chromosome("4")
                        .chromosomeBand("p1.12")
                        .gene("NF1")
                        .transcript("ENST00000358273")
                        .isCanonical(true)
                        .build())
                .build();
    }

    public static VirusInterpreterData createVirusInterpreterData()
    {
        VirusInterpreterEntry virus1 = ImmutableVirusInterpreterEntry.builder()
                .name("Human papillomavirus 16")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(1)
                .interpretation(VirusInterpretation.HPV)
                .reported(true)
                .driverLikelihood(VirusLikelihoodType.HIGH)
                .percentageCovered(0.9)
                .meanCoverage(0)
                .build();

        VirusInterpreterEntry virus2 = ImmutableVirusInterpreterEntry.builder()
                .name("Human betaherpesvirus 6B")
                .qcStatus(VirusBreakendQCStatus.NO_ABNORMALITIES)
                .integrations(0)
                .interpretation(null)     // nullable field
                .reported(false)
                .driverLikelihood(VirusLikelihoodType.LOW)
                .percentageCovered(0.4)
                .meanCoverage(0)
                .build();

        return ImmutableVirusInterpreterData.builder()
                .allViruses(List.of(virus1, virus2))
                .reportableViruses(List.of(virus1))
                .build();
    }

    public static CuppaData createCuppaData()
    {
        CuppaPrediction prediction = ImmutableCuppaPrediction.builder()
                .cancerType("Melanoma")
                .likelihood(0.996)
                .snvPairwiseClassifier(0.979)
                .genomicPositionClassifier(0.990)
                .featureClassifier(0.972)
                .altSjCohortClassifier(null)
                .expressionPairwiseClassifier(null)
                .build();

        return ImmutableCuppaData.builder()
                .mode(CuppaMode.WGS)
                .addPredictions(prediction)
                .bestPrediction(prediction)
                .simpleDups32To200B(0)
                .maxComplexSize(0)
                .telomericSGLs(0)
                .lineCount(0)
                .build();
    }

    public static LilacRecord createLilacRecord()
    {
        LilacAllele lilacAllele = ImmutableLilacAllele.builder()
                .allele("A*01:01")
                .tumorCopyNumber(6.1)
                .somaticMissense(5.0)
                .somaticNonsenseOrFrameshift(4.0)
                .somaticSplice(3.0)
                .somaticSynonymous(0.4)
                .somaticInframeIndel(1.0)
                .refFragments(210)
                .tumorFragments(1602)
                .rnaFragments(0)
                .build();
        return ImmutableLilacRecord.builder()
                .qc("PASS")
                .addAlleles(lilacAllele)
                .build();
    }

    public static ImmuneEscapeRecord createImmuneEscapeRecord()
    {
        return ImmutableImmuneEscapeRecord.builder()
                .hasHlaEscape(true)
                .hasAntigenPresentationPathwayEscape(false)
                .hasIFNGammaPathwayEscape(false)
                .hasPDL1OverexpressionEscape(false)
                .hasCD58InactivationEscape(false)
                .hasEpigeneticSETDB1Escape(false)
                .build();
    }
}