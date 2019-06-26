package com.hartwig.hmftools.patientreporter;

import java.io.IOException;
import java.util.List;
import java.util.Optional;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocation;
import com.hartwig.hmftools.common.ecrf.projections.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingChoice;
import com.hartwig.hmftools.common.purple.gender.Gender;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.purity.FittedPurityFile;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.variant.ClonalityCutoffKernel;
import com.hartwig.hmftools.common.variant.ClonalityFactory;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariantFactory;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;
import com.hartwig.hmftools.patientreporter.actionability.ClinicalTrialFactory;
import com.hartwig.hmftools.patientreporter.actionability.ReportableEvidenceItemFactory;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalysis;
import com.hartwig.hmftools.patientreporter.copynumber.CopyNumberAnalyzer;
import com.hartwig.hmftools.patientreporter.structural.SvAnalysis;
import com.hartwig.hmftools.patientreporter.structural.SvAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariant;
import com.hartwig.hmftools.patientreporter.variants.ReportableVariantAnalyzer;
import com.hartwig.hmftools.patientreporter.variants.germline.BachelorFile;
import com.hartwig.hmftools.patientreporter.variants.germline.FilterGermlineVariants;
import com.hartwig.hmftools.patientreporter.variants.germline.GermlineVariant;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalysis;
import com.hartwig.hmftools.patientreporter.variants.somatic.SomaticVariantAnalyzer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

class PatientReporter {
    private static final Logger LOGGER = LogManager.getLogger(PatientReporter.class);

    @NotNull
    private final BaseReportData baseReportData;
    @NotNull
    private final SequencedReportData sequencedReportData;
    @NotNull
    private final SvAnalyzer svAnalyzer;

    PatientReporter(@NotNull final BaseReportData baseReportData, @NotNull final SequencedReportData sequencedReportData,
            @NotNull final SvAnalyzer svAnalyzer) {
        this.baseReportData = baseReportData;
        this.sequencedReportData = sequencedReportData;
        this.svAnalyzer = svAnalyzer;
    }

    @NotNull
    public AnalysedPatientReport run(@NotNull String tumorSample, @NotNull String refSample, @NotNull String purplePurityTsv,
            @NotNull String purpleGeneCnvTsv, @NotNull String somaticVariantVcf, @Nullable String bachelorCsv,
            @NotNull String chordPredictionFile, @NotNull String circosFile, @Nullable String comments) throws IOException {
        final PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findPatientTumorLocationForSample(baseReportData.patientTumorLocations(), tumorSample);

        final CopyNumberAnalysis copyNumberAnalysis = analyzeCopyNumbers(purplePurityTsv, purpleGeneCnvTsv, patientTumorLocation);
        final SomaticVariantAnalysis somaticVariantAnalysis = analyzeSomaticVariants(tumorSample,
                somaticVariantVcf,
                copyNumberAnalysis.gender(),
                copyNumberAnalysis.purity(),
                patientTumorLocation);
        final List<GermlineVariant> germlineVariantsToReport =
                analyzeGermlineVariants(tumorSample, bachelorCsv, copyNumberAnalysis, somaticVariantAnalysis);

        final List<ReportableVariant> reportableVariants =
                ReportableVariantAnalyzer.mergeSomaticAndGermlineVariants(somaticVariantAnalysis.variantsToReport(),
                        somaticVariantAnalysis.driverCatalog(),
                        sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                        sequencedReportData.panelGeneModel().drupActionableGenes(),
                        germlineVariantsToReport,
                        sequencedReportData.germlineReportingModel(),
                        baseReportData.limsModel().germlineReportingChoice(tumorSample));

        final SvAnalysis svAnalysis = analyzeStructuralVariants(copyNumberAnalysis, patientTumorLocation);
        final ChordAnalysis chordAnalysis = analyzeChord(chordPredictionFile);

        final String clinicalSummary = sequencedReportData.summaryModel().findSummaryForSample(tumorSample);

        LOGGER.info("Printing analysis results for {}:", tumorSample);
        LOGGER.info(" Clinical summary present: " + (!clinicalSummary.isEmpty() ? "yes" : "no"));
        LOGGER.info(" Somatic variants to report : " + somaticVariantAnalysis.variantsToReport().size());
        LOGGER.info(" Germline variants to report: " + germlineVariantsToReport.size());
        LOGGER.info("  Total number of reportable variants: " + reportableVariants.size());
        LOGGER.info(" Microsatellite Indels per Mb: " + somaticVariantAnalysis.microsatelliteIndelsPerMb());
        LOGGER.info(" Tumor mutational load: " + somaticVariantAnalysis.tumorMutationalLoad());
        LOGGER.info(" Tumor mutational burden: " + somaticVariantAnalysis.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: " + chordAnalysis.hrdValue());
        LOGGER.info(" Copy number events to report: " + copyNumberAnalysis.reportableGeneCopyNumbers().size());
        LOGGER.info(" Gene fusions to report : " + svAnalysis.reportableFusions().size());
        LOGGER.info(" Gene disruptions to report : " + svAnalysis.reportableDisruptions().size());
        LOGGER.info("Printing actionability results (including off-label trials):");
        LOGGER.info(" Evidence items found based on variants: " + somaticVariantAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on copy numbers: " + copyNumberAnalysis.evidenceItems().size());
        LOGGER.info(" Evidence items found based on fusions: " + svAnalysis.evidenceItems().size());

        final List<EvidenceItem> allEvidenceItems = Lists.newArrayList();
        allEvidenceItems.addAll(somaticVariantAnalysis.evidenceItems());
        allEvidenceItems.addAll(copyNumberAnalysis.evidenceItems());
        allEvidenceItems.addAll(svAnalysis.evidenceItems());

        final SampleReport sampleReport = SampleReportFactory.fromLimsAndHospitalModel(tumorSample,
                refSample,
                baseReportData.limsModel(),
                baseReportData.hospitalModel(),
                patientTumorLocation);

        final List<EvidenceItem> nonTrials = ReportableEvidenceItemFactory.extractNonTrials(allEvidenceItems);

        return ImmutableAnalysedPatientReport.of(sampleReport,
                copyNumberAnalysis.hasReliablePurityFit(),
                copyNumberAnalysis.purity(),
                copyNumberAnalysis.ploidy(),
                clinicalSummary,
                nonTrials.stream().filter(EvidenceItem::isOnLabel).collect(Collectors.toList()),
                ClinicalTrialFactory.extractOnLabelTrials(allEvidenceItems),
                nonTrials.stream().filter(item -> !item.isOnLabel()).collect(Collectors.toList()),
                reportableVariants,
                somaticVariantAnalysis.microsatelliteIndelsPerMb(),
                somaticVariantAnalysis.tumorMutationalLoad(),
                somaticVariantAnalysis.tumorMutationalBurden(),
                chordAnalysis,
                copyNumberAnalysis.reportableGeneCopyNumbers(),
                svAnalysis.reportableFusions(),
                svAnalysis.reportableDisruptions(),
                circosFile,
                Optional.ofNullable(comments),
                baseReportData.signaturePath(),
                baseReportData.logoRVAPath(),
                baseReportData.logoCompanyPath());
    }

    @NotNull
    private CopyNumberAnalysis analyzeCopyNumbers(@NotNull String purplePurityTsv, @NotNull String purpleGeneCnvTsv,
            @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        LOGGER.info("Loading purple purity data from {}", purplePurityTsv);
        final PurityContext purityContext = FittedPurityFile.read(purplePurityTsv);
        LOGGER.info(" Purple purity " + purityContext.bestFit().purity());
        LOGGER.info(" Purple average tumor ploidy: " + purityContext.bestFit().ploidy());
        LOGGER.info(" Purple status " + purityContext.status());

        LOGGER.info("Loading purple gene copynumber data from {}", purpleGeneCnvTsv);
        final List<GeneCopyNumber> exomeGeneCopyNumbers = GeneCopyNumberFile.read(purpleGeneCnvTsv);
        LOGGER.info(" Loaded {} gene copy numbers", exomeGeneCopyNumbers.size());

        LOGGER.info("Analyzing purple copy numbers");
        return CopyNumberAnalyzer.run(purityContext,
                exomeGeneCopyNumbers,
                sequencedReportData.panelGeneModel(),
                sequencedReportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private SomaticVariantAnalysis analyzeSomaticVariants(@NotNull String sample, @NotNull String somaticVariantVcf, @NotNull Gender gender,
            double purity, @Nullable PatientTumorLocation patientTumorLocation) throws IOException {
        LOGGER.info("Loading somatic variants from {}", somaticVariantVcf);
        final List<SomaticVariant> variants = SomaticVariantFactory.filteredInstanceWithEnrichment(new PassingVariantFilter(),
                sequencedReportData.somaticVariantEnrichment()).fromVCFFile(sample, somaticVariantVcf);
        LOGGER.info(" {} PASS somatic variants loaded", variants.size());

        LOGGER.info("Enriching somatic variants");
        final List<EnrichedSomaticVariant> enrichedSomaticVariants =
                enrich(variants, gender, purity, sequencedReportData.highConfidenceRegions(), sequencedReportData.refGenomeFastaFile());

        LOGGER.info("Analyzing somatic variants");
        return SomaticVariantAnalyzer.run(enrichedSomaticVariants,
                sequencedReportData.panelGeneModel().somaticVariantGenes(),
                sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                sequencedReportData.actionabilityAnalyzer(),
                patientTumorLocation);
    }

    @NotNull
    private static List<EnrichedSomaticVariant> enrich(@NotNull List<SomaticVariant> variants, @NotNull Gender gender, double purity,
            @NotNull Multimap<String, GenomeRegion> highConfidenceRegions, @NotNull IndexedFastaSequenceFile refGenomeFastaFile) {
        final double clonalPloidy = ClonalityCutoffKernel.clonalCutoff(variants);
        final ClonalityFactory clonalityFactory = new ClonalityFactory(gender, purity, clonalPloidy);

        final EnrichedSomaticVariantFactory enrichedSomaticFactory =
                new EnrichedSomaticVariantFactory(highConfidenceRegions, refGenomeFastaFile, clonalityFactory);

        return enrichedSomaticFactory.enrich(variants);
    }

    @NotNull
    private List<GermlineVariant> analyzeGermlineVariants(@NotNull String sample, @Nullable String bachelorCsv,
            @NotNull CopyNumberAnalysis copyNumberAnalysis, @NotNull SomaticVariantAnalysis somaticVariantAnalysis) throws IOException {
        if (bachelorCsv == null) {
            LOGGER.info("Skipping germline analysis - no bachelor CSV passed. Presumably no pathogenic germline variants found.");
            return Lists.newArrayList();
        }

        LOGGER.info("Loading germline variants from {}", bachelorCsv);
        List<GermlineVariant> variants =
                BachelorFile.loadBachelorFile(bachelorCsv).stream().filter(GermlineVariant::passFilter).collect(Collectors.toList());
        LOGGER.info(" {} PASS germline variants loaded", variants.size());

        LimsGermlineReportingChoice germlineChoice = baseReportData.limsModel().germlineReportingChoice(sample);
        if (germlineChoice == LimsGermlineReportingChoice.UNKNOWN) {
            LOGGER.info(" No germline reporting choice known. No germline variants will be reported!");
            return Lists.newArrayList();
        } else {
            return FilterGermlineVariants.filterGermlineVariantsForReporting(variants,
                    sequencedReportData.germlineReportingModel(),
                    sequencedReportData.panelGeneModel().geneDriverCategoryMap(),
                    copyNumberAnalysis.exomeGeneCopyNumbers(),
                    somaticVariantAnalysis.variantsToReport());
        }
    }

    @NotNull
    private SvAnalysis analyzeStructuralVariants(@NotNull CopyNumberAnalysis copyNumberAnalysis,
            @Nullable PatientTumorLocation patientTumorLocation) {
        return svAnalyzer.run(copyNumberAnalysis.exomeGeneCopyNumbers(), sequencedReportData.actionabilityAnalyzer(), patientTumorLocation);
    }

    @NotNull
    private static ChordAnalysis analyzeChord(@NotNull String chordPredictionFile) throws IOException {
        LOGGER.info("Loading CHORD analysis from {}", chordPredictionFile);
        return ChordFileReader.read(chordPredictionFile);
    }
}
