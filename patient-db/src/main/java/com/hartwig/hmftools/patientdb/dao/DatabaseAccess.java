package com.hartwig.hmftools.patientdb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.amber.AmberAnonymous;
import com.hartwig.hmftools.common.amber.AmberMapping;
import com.hartwig.hmftools.common.amber.AmberPatient;
import com.hartwig.hmftools.common.amber.AmberSample;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsMutationalLoad;
import com.hartwig.hmftools.common.drivercatalog.dnds.DndsVariant;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.hla.HlaType;
import com.hartwig.hmftools.common.hla.HlaTypeDetails;
import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.peach.PeachCalls;
import com.hartwig.hmftools.common.peach.PeachGenotype;
import com.hartwig.hmftools.common.protect.ProtectEvidence;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.purity.SamplePurity;
import com.hartwig.hmftools.common.purple.PurpleQC;
import com.hartwig.hmftools.common.sigs.SignatureAllocation;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxBreakend;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxFusion;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvAnnotation;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertion;
import com.hartwig.hmftools.patientdb.clinical.datamodel.Patient;
import com.hartwig.hmftools.patientdb.clinical.datamodel.SampleData;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.conf.MappedSchema;
import org.jooq.conf.RenderMapping;
import org.jooq.conf.Settings;
import org.jooq.impl.DSL;

import htsjdk.variant.variantcontext.VariantContext;

public class DatabaseAccess implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(DatabaseAccess.class);
    private static final String DEV_CATALOG = "hmfpatients_test";

    public static final double MIN_SAMPLE_PURITY = 0.195;

    public static final String DB_USER = "db_user";
    public static final String DB_PASS = "db_pass";
    public static final String DB_URL = "db_url";

    public static final String DB_DEFAULT_ARGS = "?serverTimezone=UTC&useSSL=false";

    @NotNull
    private final DSLContext context;
    @NotNull
    private final EcrfDAO ecrfDAO;
    @NotNull
    private final ClinicalDAO clinicalDAO;
    @NotNull
    private final ValidationFindingDAO validationFindingsDAO;
    @NotNull
    private final RNADAO rnaDAO;
    @NotNull
    private final PurityDAO purityDAO;
    @NotNull
    private final AmberDAO amberDAO;
    @NotNull
    private final MetricDAO metricDAO;
    @NotNull
    private final PeachDAO pgxDAO;
    @NotNull
    private final CopyNumberDAO copyNumberDAO;
    @NotNull
    private final DriverCatalogDAO driverCatalogDAO;
    @NotNull
    private final GeneCopyNumberDAO geneCopyNumberDAO;
    @NotNull
    private final SomaticVariantDAO somaticVariantDAO;
    @NotNull
    private final StructuralVariantDAO structuralVariantDAO;
    @NotNull
    private final StructuralVariantClusterDAO structuralVariantClusterDAO;
    @NotNull
    private final StructuralVariantFusionDAO structuralVariantFusionDAO;
    @NotNull
    private final SignatureDAO signatureDAO;
    @NotNull
    private final CanonicalTranscriptDAO canonicalTranscriptDAO;
    @NotNull
    private final ChordDAO chordDAO;
    @NotNull
    private final ProtectDAO protectDAO;
    @NotNull
    private final DriverGenePanelDAO driverGenePanelDAO;
    @NotNull
    private final GermlineVariantDAO germlineVariantDAO;
    @NotNull
    private final HlaTypeDAO hlaTypeDAO;

    public DatabaseAccess(@NotNull final String userName, @NotNull final String password, @NotNull final String url) throws SQLException {
        // Disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        Connection conn = DriverManager.getConnection(url, userName, password);
        String catalog = conn.getCatalog();
        LOGGER.debug("Connecting to database {}", catalog);
        this.context = DSL.using(conn, SQLDialect.MYSQL, settings(catalog));

        ecrfDAO = new EcrfDAO(context);
        clinicalDAO = new ClinicalDAO(context);
        protectDAO = new ProtectDAO(context);
        validationFindingsDAO = new ValidationFindingDAO(context);
        rnaDAO = new RNADAO(context);
        purityDAO = new PurityDAO(context);
        amberDAO = new AmberDAO(context);
        metricDAO = new MetricDAO(context);
        pgxDAO = new PeachDAO(context);
        copyNumberDAO = new CopyNumberDAO(context);
        geneCopyNumberDAO = new GeneCopyNumberDAO(context);
        somaticVariantDAO = new SomaticVariantDAO(context);
        structuralVariantDAO = new StructuralVariantDAO(context);
        structuralVariantClusterDAO = new StructuralVariantClusterDAO(context);
        structuralVariantFusionDAO = new StructuralVariantFusionDAO(context);
        signatureDAO = new SignatureDAO(context);
        canonicalTranscriptDAO = new CanonicalTranscriptDAO(context);
        driverCatalogDAO = new DriverCatalogDAO(context);
        chordDAO = new ChordDAO(context);
        driverGenePanelDAO = new DriverGenePanelDAO(context);
        germlineVariantDAO = new GermlineVariantDAO(context);
        hlaTypeDAO = new HlaTypeDAO(context);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options) {
        addDatabaseCmdLineArgs(options, false);
    }

    public static void addDatabaseCmdLineArgs(@NotNull Options options, boolean isRequired) {
        options.addOption(Option.builder(DB_USER).desc("Database username").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_PASS).desc("Database password").hasArg(true).required(isRequired).build());
        options.addOption(Option.builder(DB_URL).desc("Database url").hasArg(true).required(isRequired).build());
    }

    public static boolean hasDatabaseConfig(@NotNull CommandLine cmd) {
        return cmd.hasOption(DB_URL) && cmd.hasOption(DB_USER) && cmd.hasOption(DB_PASS);
    }

    @NotNull
    public static DatabaseAccess databaseAccess(@NotNull CommandLine cmd) throws SQLException {
        return databaseAccess(cmd, false);
    }

    @NotNull
    public static DatabaseAccess databaseAccess(@NotNull CommandLine cmd, boolean applyDefaultArgs) throws SQLException {
        String userName = cmd.getOptionValue(DB_USER);
        String password = cmd.getOptionValue(DB_PASS);
        String databaseUrl = cmd.getOptionValue(DB_URL);
        String jdbcUrl = "jdbc:" + databaseUrl;

        if (applyDefaultArgs && !jdbcUrl.contains("serverTimezone") && !jdbcUrl.contains("useSSL")) {
            jdbcUrl += DB_DEFAULT_ARGS;
        }

        return new DatabaseAccess(userName, password, jdbcUrl);
    }

    @Nullable
    public static DatabaseAccess createDatabaseAccess(@NotNull CommandLine cmd) {
        if (!hasDatabaseConfig(cmd)) {
            return null;
        }

        try {
            return databaseAccess(cmd, true);
        } catch (SQLException e) {
            LOGGER.error("DB connection failed: {}", e.toString());
            return null;
        }
    }

    @NotNull
    public DSLContext context() {
        return context;
    }

    @Override
    public void close() {
        context.close();
    }

    @Nullable
    private static Settings settings(@NotNull String catalog) {
        return !catalog.equals(DEV_CATALOG)
                ? new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)))
                : null;
    }

    @NotNull
    public BufferedWriter<VariantContext> germlineVariantWriter(String tumorSample, String referenceSample, String rnaSample) {
        return germlineVariantDAO.writer(tumorSample, referenceSample, rnaSample);
    }

    @NotNull
    public List<String> readPurpleSampleList() {
        return purityDAO.getSampleIds();
    }

    @NotNull
    public List<String> readPurpleSampleListPassingQC(double minPurity) {
        return purityDAO.getSamplesPassingQC(minPurity);
    }

    @NotNull
    public List<SamplePurity> readSamplePurityPassingQC(double minPurity) {
        return purityDAO.readPassingQC(minPurity);
    }

    @NotNull
    public List<SamplePurity> readSamplePurityPassingQC(double minPurity, @NotNull String primaryTumorLocation) {
        return purityDAO.readPassingQC(minPurity, primaryTumorLocation);
    }

    @Nullable
    public PurityContext readPurityContext(@NotNull String sampleId) {
        return purityDAO.readPurityContext(sampleId);
    }

    @NotNull
    public List<PurpleCopyNumber> readCopynumbers(@NotNull String sample) {
        return copyNumberDAO.read(sample);
    }

    @NotNull
    public List<GeneCopyNumber> readGeneCopynumbers(@NotNull String sample, @NotNull List<String> genes) {
        return geneCopyNumberDAO.read(sample, genes);
    }

    @NotNull
    public List<String> readSomaticVariantSampleList() {
        return somaticVariantDAO.getSamplesList();
    }

    @NotNull
    public List<DndsVariant> readDndsVariants(int maxRepeatCount, @NotNull String sample) {
        return somaticVariantDAO.readDndsVariants(maxRepeatCount, sample);
    }

    @NotNull
    public DndsMutationalLoad readDndsMutationLoad(@NotNull String sample) {
        return somaticVariantDAO.readDndsLoad(sample);
    }

    @NotNull
    public List<SomaticVariant> readSomaticVariants(@NotNull String sample, VariantType type) {
        return somaticVariantDAO.read(sample, type);
    }

    @NotNull
    public List<String> readStructuralVariantSampleList(@NotNull String sampleSearch) {
        return structuralVariantDAO.getSamplesList(sampleSearch);
    }

    @NotNull
    public List<StructuralVariantData> readStructuralVariantData(@NotNull String sample) {
        return structuralVariantDAO.read(sample);
    }

    @NotNull
    public List<LinxCluster> readClusters(@NotNull String sample) {
        return structuralVariantClusterDAO.readClusters(sample);
    }

    @NotNull
    public List<LinxSvAnnotation> readSvAnnotations(@NotNull String sample) {
        return structuralVariantClusterDAO.readAnnotations(sample);
    }

    @NotNull
    public List<DriverCatalog> readDriverCatalog(@NotNull String sample) {
        return driverCatalogDAO.readDriverData(sample);
    }

    @NotNull
    public List<LinxDriver> readSvDriver(@NotNull String sample) {
        return structuralVariantClusterDAO.readSvDrivers(sample);
    }

    @NotNull
    public List<SignatureAllocation> readSignatureAllocations(@NotNull String sample) {
        return signatureDAO.readAllocations(sample);
    }

    @NotNull
    public List<LinxFusion> readFusions(@NotNull String sample) {
        return structuralVariantFusionDAO.readFusions(sample);
    }

    @NotNull
    public List<LinxBreakend> readBreakends(@NotNull String sample) {
        return structuralVariantFusionDAO.readBreakends(sample);
    }

    @NotNull
    public List<LinxViralInsertion> readViralInsertions(@NotNull String sample) {
        return structuralVariantFusionDAO.readViralInsertions(sample);
    }

    public void writeCanonicalTranscripts(@NotNull String assembly, @NotNull List<CanonicalTranscript> transcripts) {
        canonicalTranscriptDAO.write(assembly, transcripts);
    }

    public void writePurity(@NotNull String sampleId, @NotNull PurityContext context, @NotNull PurpleQC checks) {
        purityDAO.write(sampleId, context, checks);
    }

    public void writeBestFitPerPurity(@NotNull String sampleId, @NotNull List<FittedPurity> bestFitPerPurity) {
        purityDAO.write(sampleId, bestFitPerPurity);
    }

    public void writeCopynumbers(@NotNull String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeCopyNumber(sample, copyNumbers);
    }

    public void writeAmberMapping(@NotNull String sample, List<AmberMapping> mapping) {
        amberDAO.writeMapping(sample, mapping);
    }

    public void writeAmberPatients(List<AmberPatient> mapping) {
        amberDAO.writePatients(mapping);
    }

    public void writeAmberAnonymous(List<AmberAnonymous> mapping) {
        amberDAO.writeAnonymous(mapping);
    }

    public void writeAmberSample(@NotNull AmberSample identity) {
        amberDAO.writeIdentity(identity);
    }

    @NotNull
    public List<AmberSample> readAmberSamples() {
        return amberDAO.readSamples();
    }

    @NotNull
    public List<AmberAnonymous> readAmberAnonymous() {
        return amberDAO.readAnonymous();
    }

    @NotNull
    public List<AmberPatient> readAmberPatients() {
        return amberDAO.readPatients();
    }

    public void truncateAmberPatients() {
        amberDAO.truncatePatients();
    }

    public void truncateAmberMappings() {
        amberDAO.truncateMappings();
    }

    @NotNull
    public BufferedWriter<SomaticVariant> somaticVariantWriter(@NotNull final String sampleId) {
        return somaticVariantDAO.writer(sampleId);
    }

    public void writeStructuralVariants(@NotNull String sampleId, @NotNull List<StructuralVariantData> variants) {
        structuralVariantDAO.write(sampleId, variants);
    }

    public void writeGermlineCopynumbers(@NotNull String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeGermlineCopyNumber(sample, copyNumbers);
    }

    public void writeSvClusters(@NotNull String sample, @NotNull List<LinxCluster> clusters) {
        structuralVariantClusterDAO.writeClusters(sample, clusters);
    }

    public void writeSvLinxData(@NotNull String sample, @NotNull List<LinxSvAnnotation> svData) {
        structuralVariantClusterDAO.writeSvData(sample, svData);
    }

    public void writeSvLinks(@NotNull String sample, @NotNull List<LinxLink> links) {
        structuralVariantClusterDAO.writeLinks(sample, links);
    }

    public void writeSvDrivers(@NotNull String sample, @NotNull List<LinxDriver> drivers) {
        structuralVariantClusterDAO.writeDrivers(sample, drivers);
    }

    public void writeSvViralInserts(@NotNull String sample, @NotNull List<LinxViralInsertion> inserts) {
        structuralVariantClusterDAO.writeViralInserts(sample, inserts);
    }

    public void writeSignatures(@NotNull String sample, @NotNull List<SignatureAllocation> sigAllocations) {
        signatureDAO.write(sample, sigAllocations);
    }

    public void writeGeneCopynumberRegions(@NotNull String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        geneCopyNumberDAO.writeCopyNumber(sample, geneCopyNumbers);
    }

    public void writeProtectDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> driverCatalog) {
        driverCatalogDAO.writeAll(sample, driverCatalog);
    }

    public void writeLinxDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> somaticCatalog) {
        driverCatalogDAO.writeLinx(sample, somaticCatalog);
    }

    public void writePurpleDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> somaticCatalog,
            @NotNull List<DriverCatalog> germlineCatalog) {
        driverCatalogDAO.writePurple(sample, somaticCatalog, germlineCatalog);
    }

    public void writeGermlineDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> driverCatalog) {
        driverCatalogDAO.writeGermline(sample, driverCatalog);
    }

    public void writeMetrics(@NotNull String sample, @NotNull WGSMetricWithQC metrics) {
        metricDAO.writeMetrics(sample, metrics);
    }

    public void writePeach(@NotNull String sample, @NotNull List<PeachGenotype> peachGenotypes, @NotNull List<PeachCalls> peachCalls) {
        pgxDAO.writePeach(sample, peachGenotypes, peachCalls);
    }

    public void writeProtectEvidence(@NotNull String sample, @NotNull List<ProtectEvidence> evidence) {
        protectDAO.write(sample, evidence);
    }

    public void writeChord(@NotNull String sample, @NotNull ChordAnalysis chordAnalysis) {
        chordDAO.writeChord(sample, chordAnalysis);
    }

    public void writeRNA(@NotNull Set<String> samples) {
        rnaDAO.write(samples);
    }

    public void clearCpctEcrf() {
        ecrfDAO.clearCpct();
    }

    public void clearDrupEcrf() {
        ecrfDAO.clearDrup();
    }

    public void clearClinicalTables() {
        validationFindingsDAO.clear();
        clinicalDAO.clear();
    }

    public void writeHla(@NotNull final String sample, @NotNull final HlaType type, @NotNull final List<HlaTypeDetails> details) {
        hlaTypeDAO.writeType(sample, type);
        hlaTypeDAO.writeTypeDetails(sample, details);
    }

    public void writeFullClinicalData(@NotNull Patient patient, boolean blacklisted) {
        clinicalDAO.writeFullClinicalData(patient, blacklisted);
    }

    public void writeSampleClinicalData(@NotNull String patientIdentifier, boolean blacklisted, @NotNull List<SampleData> samples) {
        clinicalDAO.writeSampleClinicalData(patientIdentifier, blacklisted, samples);
    }

    public void writeGenePanel(DriverGenePanel panel) {
        driverGenePanelDAO.writePanel(panel);
    }

    public void writeDrupEcrf(@NotNull EcrfModel model, @NotNull Set<String> sequencedPatients) {
        LOGGER.info(" Writing DRUP datamodel...");
        ecrfDAO.writeDrupDatamodel(model.fields());
        LOGGER.info("  Done writing DRUP datamodel.");
        LOGGER.info(" Writing raw DRUP patient data...");
        model.patients().forEach(patient -> ecrfDAO.writeDrupPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("  Done writing raw DRUP patient data.");
    }

    public void writeCpctEcrf(@NotNull EcrfModel model, @NotNull Set<String> sequencedPatients) {
        LOGGER.info(" Writing CPCT datamodel...");
        ecrfDAO.writeCpctDatamodel(model.fields());
        LOGGER.info("  Done writing CPCT datamodel.");
        LOGGER.info(" Writing raw CPCT patient data...");
        model.patients().forEach(patient -> ecrfDAO.writeCpctPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("  Done writing raw CPCT patient data.");
    }

    public void writeValidationFindings(@NotNull List<ValidationFinding> findings) {
        validationFindingsDAO.write(findings);
    }

    public void deleteAllDataForSample(@NotNull String sample) {
        LOGGER.info("Deleting metric data for sample: {}", sample);
        metricDAO.deleteMetricForSample(sample);

        LOGGER.info("Deleting CHORD data for sample: {}", sample);
        chordDAO.deleteChordForSample(sample);

        LOGGER.info("Deleting AMBER data for sample: {}", sample);
        amberDAO.deleteAmberRecordsForSample(sample);

        LOGGER.info("Deleting purity data for sample: {}", sample);
        purityDAO.deletePurityForSample(sample);

        LOGGER.info("Deleting copy number data for sample: {}", sample);
        copyNumberDAO.deleteCopyNumberForSample(sample);

        LOGGER.info("Deleting gene copy numbers for sample: {}", sample);
        geneCopyNumberDAO.deleteGeneCopyNumberForSample(sample);

        LOGGER.info("Deleting somatic variants for sample: {}", sample);
        somaticVariantDAO.deleteSomaticVariantForSample(sample);

        LOGGER.info("Deleting germline variant data for sample: {}", sample);
        context.delete(Tables.GERMLINEVARIANT).where(Tables.GERMLINEVARIANT.SAMPLEID.eq(sample)).execute();
        germlineVariantDAO.deleteGermlineVariantsForSample(sample);

        LOGGER.info("Deleting structural variant annotation data for sample: {}", sample);
        structuralVariantFusionDAO.deleteAnnotationsForSample(sample);

        LOGGER.info("Deleting structural variant cluster data for sample: {}", sample);
        structuralVariantClusterDAO.deleteClusterDataForSample(sample);

        LOGGER.info("Deleting structural variants for sample: {}", sample);
        structuralVariantDAO.deleteStructuralVariantsForSample(sample);

        LOGGER.info("Deleting signatures for sample: {}", sample);
        signatureDAO.deleteSignatureDataForSample(sample);

        LOGGER.info("Deleting PROTECT data for sample: {}", sample);
        protectDAO.deleteEvidenceForSample(sample);

        LOGGER.info("Deleting driver catalog for sample: {}", sample);
        driverCatalogDAO.deleteForSample(sample);

        LOGGER.info("Deleting pgx data for sample: {}", sample);
        pgxDAO.deletePgxForSample(sample);

        LOGGER.info("Deleting hla data for sample: {}", sample);
        hlaTypeDAO.deleteHlaFprSample(sample);

        LOGGER.info("All data for sample '{}' has been deleted", sample);
    }
}

