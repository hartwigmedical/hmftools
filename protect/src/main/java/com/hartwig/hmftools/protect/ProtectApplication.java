package com.hartwig.hmftools.protect;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.EvidenceItemFile;
import com.hartwig.hmftools.common.clinical.PatientTumorLocation;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFile;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationFunctions;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationV2;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationV2File;
import com.hartwig.hmftools.common.clinical.PatientTumorLocationV2Functions;
import com.hartwig.hmftools.common.doid.DiseaseOntology;
import com.hartwig.hmftools.common.doid.DoidParents;
import com.hartwig.hmftools.common.lims.LimsGermlineReportingLevel;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItem;
import com.hartwig.hmftools.common.protect.ProtectEvidenceItemFile;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;
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

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class ProtectApplication implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(ProtectApplication.class);
    private static final String PAN_CANCER_DOID = "162";

    public static void main(@NotNull String[] args) throws IOException {
        Options options = ProtectConfig.createOptions();
        DatabaseAccess.addDatabaseCmdLineArgs(options);

        try (final ProtectApplication application = new ProtectApplication(options, args)) {
            application.run();
            application.runOld();
        } catch (ParseException exception) {
            LOGGER.warn(exception);
            new HelpFormatter().printHelp("PROTECT", options);
            System.exit(1);
        } catch (SQLException exception) {
            LOGGER.warn(exception);
            System.exit(1);
        }
    }

    private final DatabaseAccess dbAccess;
    private final ProtectConfig protectConfig;

    public ProtectApplication(final Options options, final String... args) throws ParseException, SQLException, IOException {
        final CommandLine cmd = new DefaultParser().parse(options, args);
        this.dbAccess = DatabaseAccess.databaseAccess(cmd);
        this.protectConfig = ProtectConfig.createConfig(cmd);
    }

    public void run() throws IOException {
        final List<ProtectEvidenceItem> evidence = protectEvidence(protectConfig);

        LOGGER.info("Writing {} records to database", evidence.size());
        dbAccess.writeProtectEvidence(protectConfig.tumorSampleId(), evidence);

        final String filename = ProtectEvidenceItemFile.generateFilename(protectConfig.outputDir(), protectConfig.tumorSampleId());
        LOGGER.info("Writing {} records to file: {}", evidence.size(), filename);
        ProtectEvidenceItemFile.write(filename, evidence);
    }

    public void runOld() throws IOException {
        String tumorSampleId = protectConfig.tumorSampleId();
        LOGGER.info("Running PROTECT for {}", tumorSampleId);

        PatientTumorLocation patientTumorLocation = loadPatientTumorLocation(protectConfig.tumorLocationTsvV1(), tumorSampleId);
        LOGGER.info("Creating deprecated actionability analyzer from {}", protectConfig.deprecatedActionabilityDir());
        ActionabilityAnalyzer actionabilityAnalyzer = ActionabilityAnalyzer.fromKnowledgebase(protectConfig.deprecatedActionabilityDir());
        LOGGER.info("Creating germline reporting model from {}", protectConfig.germlineGenesCsv());
        GermlineReportingModel germlineReportingModel = GermlineReportingFile.buildFromCsv(protectConfig.germlineGenesCsv());

        GenomicAnalyzer analyzer = new GenomicAnalyzer(actionabilityAnalyzer, germlineReportingModel);
        GenomicAnalysis analysis = analyzer.run(tumorSampleId,
                patientTumorLocation,
                LimsGermlineReportingLevel.REPORT_WITHOUT_NOTIFICATION,
                true,
                protectConfig.purplePurityTsv(),
                protectConfig.purpleQcFile(),
                protectConfig.purpleDriverCatalogTsv(),
                protectConfig.purpleSomaticVariantVcf(),
                protectConfig.bachelorTsv(),
                protectConfig.linxFusionTsv(),
                protectConfig.linxBreakendTsv(),
                protectConfig.linxViralInsertionTsv(),
                protectConfig.linxDriversTsv(),
                protectConfig.chordPredictionTxt());

        //        printResults(tumorSampleId, analysis);

        EvidenceItemFile.write(protectConfig.outputDir() + File.separator + protectConfig.tumorSampleId() + ".old.offLabel.tsv",
                analysis.offLabelEvidence());
        EvidenceItemFile.write(protectConfig.outputDir() + File.separator + protectConfig.tumorSampleId() + ".old.onLabel.tsv",
                analysis.tumorSpecificEvidence());
    }

    private static Set<String> doids(@NotNull ProtectConfig config) throws IOException {
        final Set<String> result = Sets.newHashSet();
        LOGGER.info("Loading DOID file from {}", config.doidJsonFile());
        final DoidParents doidParent = new DoidParents(DiseaseOntology.readDoidJsonFile(config.doidJsonFile()).doidEdges());

        LOGGER.info("Loading patient tumor locations from {}", config.tumorLocationTsvV2());
        final List<PatientTumorLocationV2> tumorLocations = PatientTumorLocationV2File.read(config.tumorLocationTsvV2());
        @Nullable
        PatientTumorLocationV2 sampleTumorLocation =
                PatientTumorLocationV2Functions.findTumorLocationForSample(tumorLocations, config.tumorSampleId());
        if (sampleTumorLocation == null) {
            result.add(PAN_CANCER_DOID);
            return result;
        }

        final List<String> initialDoids = sampleTumorLocation.doids();
        for (String initialDoid : initialDoids) {
            result.add(initialDoid);
            result.addAll(doidParent.parents(initialDoid));
        }

        LOGGER.info(" Resolved doids: {}", String.join(";", result));
        return result;
    }

    @NotNull
    private static List<ProtectEvidenceItem> protectEvidence(ProtectConfig config) throws IOException {
        final Set<String> doids = doids(config);

        // Serve Data
        final String serveActionabilityDir = config.serveActionabilityDir();
        final ActionableEvents actionableEvents = ActionableEventsLoader.readFromDir(serveActionabilityDir);
        final VariantEvidence variantEvidenceFactory = new VariantEvidence(actionableEvents.hotspots(), actionableEvents.ranges());
        final CopyNumberEvidence copyNumberEvidenceFactory = new CopyNumberEvidence(actionableEvents.genes());
        final FusionEvidence fusionEvidenceFactory = new FusionEvidence(actionableEvents.genes(), actionableEvents.fusions());

        // External Data
        final LinxData linxData = LinxDataLoader.load(config);
        final PurpleData purpleData = PurpleDataLoader.load(config);
        final BachelorData bachelorData = BachelorDataLoader.load(config.bachelorTsv(), purpleData, linxData);

        // Evidence
        final List<ProtectEvidenceItem> variantEvidence =
                variantEvidenceFactory.evidence(doids, bachelorData.germlineVariants(), purpleData.somaticVariants());
        final List<ProtectEvidenceItem> copyNumberEvidence = copyNumberEvidenceFactory.evidence(doids, purpleData.copyNumberAlterations());
        final List<ProtectEvidenceItem> fusionEvidence = fusionEvidenceFactory.evidence(doids, linxData.fusions());

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

    @Override
    public void close() {
        dbAccess.close();
        LOGGER.info("Complete");
    }
}
