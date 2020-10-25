package com.hartwig.hmftools.protect;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItemFile;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.chord.ChordFileReader;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItemFile;
import com.hartwig.hmftools.protect.bachelor.BachelorData;
import com.hartwig.hmftools.protect.bachelor.BachelorDataLoader;
import com.hartwig.hmftools.protect.evidence.CopyNumberEvidence;
import com.hartwig.hmftools.protect.evidence.FusionEvidence;
import com.hartwig.hmftools.protect.evidence.VariantEvidence;
import com.hartwig.hmftools.protect.linx.LinxData;
import com.hartwig.hmftools.protect.linx.LinxDataLoader;
import com.hartwig.hmftools.protect.purple.PurpleData;
import com.hartwig.hmftools.protect.purple.PurpleDataLoader;
import com.hartwig.hmftools.protect.variants.ReportableVariant;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingFile;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.serve.actionability.ActionableEvents;
import com.hartwig.hmftools.serve.actionability.ActionableEventsLoader;

import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ProtectApplication {

    private static final Logger LOGGER = LogManager.getLogger(ProtectApplication.class);

    public static void main(@NotNull String[] args) throws IOException {
        Options options = ProtectConfig.createOptions();

        ProtectConfig config = null;
        try {
            config = ProtectConfig.createConfig(new DefaultParser().parse(options, args));
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PROTECT", options);
            System.exit(1);
        }

        String tumorSampleId = config.tumorSampleId();
        LOGGER.info("Running PROTECT for {}", tumorSampleId);

        PatientTumorLocation patientTumorLocation = loadPatientTumorLocation(config.tumorLocationTsv(), tumorSampleId);

        LOGGER.info("Creating deprecated actionability analyzer from {}", config.deprecatedActionabilityDir());
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(config.deprecatedActionabilityDir());
        LOGGER.info("Creating germline reporting model from {}", config.germlineGenesCsv());
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(config.germlineGenesCsv());

        GenomicAnalyzer analyzer = new GenomicAnalyzer(actionabilityAnalyzer, germlineReportingModel);
        GenomicAnalysis analysis = analyzer.run(tumorSampleId,
                patientTumorLocation,
                LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                true,
                config.purplePurityTsv(),
                config.purpleQcFile(),
                config.purpleGeneCnvTsv(),
                config.purpleDriverCatalogTsv(),
                config.purpleSomaticVariantVcf(),
                config.bachelorTsv(),
                config.linxFusionTsv(),
                config.linxBreakendTsv(),
                config.linxViralInsertionTsv(),
                config.linxDriversTsv(),
                config.chordPredictionTxt());

        printResults(tumorSampleId, analysis);

        EvidenceItemFile.write(config.outputDir() + File.separator + "evidence.off.tsv", analysis.offLabelEvidence());
        EvidenceItemFile.write(config.outputDir() + File.separator + "evidence.tumor.tsv", analysis.tumorSpecificEvidence());

        List<ProtectEvidenceItem> protectEvidence = protectEvidence(config).stream().filter(ProtectEvidenceItem::reported).collect(Collectors.toList());
        ProtectEvidenceItemFile.write(config.outputDir() + File.separator + "evidence.protect.tsv", protectEvidence);

        LOGGER.info("Complete");
    }

    @NotNull
    private static List<ProtectEvidenceItem> protectEvidence(ProtectConfig config) throws IOException {

        // Serve Data
        final String serveActionabilityDir = config.serveActionabilityDir();
        final ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(serveActionabilityDir);
        final VariantEvidence variantEvidenceFactory = new VariantEvidence(actionableEvents.hotspots(), actionableEvents.ranges());
        final CopyNumberEvidence copyNumberEvidenceFactory = new CopyNumberEvidence(actionableEvents.genes());
        final FusionEvidence fusionEvidenceFactory = new FusionEvidence(actionableEvents.genes(), actionableEvents.fusions());

        // External Data
        final ChordAnalysis chordAnalysis = ChordFileReader.read(config.chordPredictionTxt());
        final LinxData linxData = LinxDataLoader.load(config);
        final PurpleData purpleData = PurpleDataLoader.load(config);
        final BachelorData bachelorData = BachelorDataLoader.load(purpleData, config.bachelorTsv(), chordAnalysis.hrStatus());

        // Evidence
        final List<ProtectEvidenceItem> variantEvidence =
                variantEvidenceFactory.evidence(Sets.newHashSet(), bachelorData.germlineVariants(), purpleData.somaticVariants());
        final List<ProtectEvidenceItem> copyNumberEvidence =
                copyNumberEvidenceFactory.evidence(Sets.newHashSet(), purpleData.copyNumberAlterations());
        final List<ProtectEvidenceItem> fusionEvidence = fusionEvidenceFactory.evidence(Sets.newHashSet(), linxData.fusions());

        final List<ProtectEvidenceItem> result = Lists.newArrayList();
        result.addAll(variantEvidence);
        result.addAll(copyNumberEvidence);
        result.addAll(fusionEvidence);
        return result;
    }

    @Nullable
    private static PatientTumorLocation loadPatientTumorLocation(@NotNull String tumorLocationTsv, @NotNull String tumorSampleId)
            throws IOException {
        List<PatientTumorLocation> patientTumorLocationList = PatientTumorLocationFile.readRecordsTSV(tumorLocationTsv);
        LOGGER.info("Loaded {} patient tumor locations from {}", patientTumorLocationList.size(), tumorLocationTsv);
        PatientTumorLocation patientTumorLocation =
                PatientTumorLocationFunctions.findTumorLocationForSample(patientTumorLocationList, tumorSampleId);
        LOGGER.info(" Resolved tumor location to '{}' for {}", patientTumorLocation, tumorSampleId);
        return patientTumorLocation;
    }

    private static void printResults(@NotNull String tumorSampleId, @NotNull GenomicAnalysis analysis) {
        List<ReportableVariant> variantsWithNotify =
                analysis.reportableVariants().stream().filter(ReportableVariant::notifyClinicalGeneticist).collect(Collectors.toList());
        LOGGER.info("Printing genomic analysis results for {}:", tumorSampleId);
        LOGGER.info(" Somatic variants to report: {}", analysis.reportableVariants().size());
        LOGGER.info("  Variants for which to notify clinical geneticist: {}", variantsWithNotify.size());
        LOGGER.info(" Microsatellite indels per Mb: {} ({})", analysis.microsatelliteIndelsPerMb(), analysis.microsatelliteStatus());
        LOGGER.info(" Tumor mutational load: {} ({})", analysis.tumorMutationalLoad(), analysis.tumorMutationalLoadStatus());
        LOGGER.info(" Tumor mutational burden: {}", analysis.tumorMutationalBurden());
        LOGGER.info(" CHORD analysis HRD prediction: {} ({})", analysis.chordHrdValue(), analysis.chordHrdStatus());
        LOGGER.info(" Number of gains and losses to report: {}", analysis.gainsAndLosses().size());
        LOGGER.info(" Gene fusions to report: {}", analysis.geneFusions().size());
        LOGGER.info(" Gene disruptions to report: {}", analysis.geneDisruptions().size());
        LOGGER.info(" Viral insertions to report: {}", analysis.viralInsertions() != null ? analysis.viralInsertions().size() : "0");

        LOGGER.info("Printing actionability results for {}", tumorSampleId);
        LOGGER.info(" Tumor-specific evidence items found: {}", analysis.tumorSpecificEvidence().size());
        LOGGER.info(" Clinical trials matched to molecular profile: {}", analysis.clinicalTrials().size());
        LOGGER.info(" Off-label evidence items found: {}", analysis.offLabelEvidence().size());

    }
}
