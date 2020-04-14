package com.hartwig.hmftools.patientdb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.actionability.EvidenceItem;
import com.hartwig.hmftools.common.amber.AmberBAF;
import com.hartwig.hmftools.common.chord.ChordAnalysis;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.genome.region.CanonicalTranscript;
import com.hartwig.hmftools.common.metrics.WGSMetricWithQC;
import com.hartwig.hmftools.common.pharmacogenetics.PGXCalls;
import com.hartwig.hmftools.common.pharmacogenetics.PGXGenotype;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.VariantType;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxCluster;
import com.hartwig.hmftools.common.variant.structural.linx.LinxDriver;
import com.hartwig.hmftools.common.variant.structural.linx.LinxLink;
import com.hartwig.hmftools.common.variant.structural.linx.LinxSvData;
import com.hartwig.hmftools.common.variant.structural.linx.LinxViralInsertFile;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.database.hmfpatients.Tables;

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

public class DatabaseAccess implements AutoCloseable {

    private static final Logger LOGGER = LogManager.getLogger(DatabaseAccess.class);
    private static final String DEV_CATALOG = "hmfpatients_test";

    public static final double MIN_SAMPLE_PURITY = 0.195;

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
    private final PgxDAO pgxDAO;
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
    private final CanonicalTranscriptDAO canonicalTranscriptDAO;
    @NotNull
    private final ChordDAO chordDAO;
    @NotNull
    private final ClinicalEvidenceDAO clinicalEvidenceDAO;

    public DatabaseAccess(@NotNull final String userName, @NotNull final String password, @NotNull final String url) throws SQLException {
        // Disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        Connection conn = DriverManager.getConnection(url, userName, password);
        String catalog = conn.getCatalog();
        LOGGER.debug("Connecting to database {}", catalog);
        this.context = DSL.using(conn, SQLDialect.MYSQL, settings(catalog));

        ecrfDAO = new EcrfDAO(context);
        clinicalDAO = new ClinicalDAO(context);
        validationFindingsDAO = new ValidationFindingDAO(context);
        rnaDAO = new RNADAO(context);
        purityDAO = new PurityDAO(context);
        amberDAO = new AmberDAO(context);
        metricDAO = new MetricDAO(context);
        pgxDAO = new PgxDAO(context);
        copyNumberDAO = new CopyNumberDAO(context);
        geneCopyNumberDAO = new GeneCopyNumberDAO(context);
        somaticVariantDAO = new SomaticVariantDAO(context);
        structuralVariantDAO = new StructuralVariantDAO(context);
        structuralVariantClusterDAO = new StructuralVariantClusterDAO(context);
        structuralVariantFusionDAO = new StructuralVariantFusionDAO(context);
        canonicalTranscriptDAO = new CanonicalTranscriptDAO(context);
        driverCatalogDAO = new DriverCatalogDAO(context);
        chordDAO = new ChordDAO(context);
        clinicalEvidenceDAO = new ClinicalEvidenceDAO(context);
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
    public List<String> readPurpleSampleList() {
        return purityDAO.getSampleIds();
    }

    @NotNull
    public List<String> readPurpleSampleListPassingQC(double minPurity) {
        return purityDAO.getSamplesPassingQC(minPurity);
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
    public List<DriverCatalog> readDriverCatalog(@NotNull String sample) {
        return driverCatalogDAO.readDriverData(sample);
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

    public void writeAmberBAF(@NotNull String sampleId, @NotNull List<AmberBAF> amber) {
        amberDAO.write(sampleId, amber);
    }

    public void writeSomaticVariants(@NotNull String sampleId, @NotNull List<SomaticVariant> variants) {
        somaticVariantDAO.write(sampleId, variants);
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

    public void writeSvLinxData(@NotNull String sample, @NotNull List<LinxSvData> svData) {
        structuralVariantClusterDAO.writeSvData(sample, svData);
    }

    public void writeSvLinks(@NotNull String sample, @NotNull List<LinxLink> links) {
        structuralVariantClusterDAO.writeLinks(sample, links);
    }

    public void writeSvDrivers(@NotNull String sample, @NotNull List<LinxDriver> drivers) {
        structuralVariantClusterDAO.writeDrivers(sample, drivers);
    }

    public void writeSvViralInserts(@NotNull String sample, @NotNull List<LinxViralInsertFile> inserts) {
        structuralVariantClusterDAO.writeViralInserts(sample, inserts);
    }

    public void writeGeneCopynumberRegions(@NotNull String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        geneCopyNumberDAO.writeCopyNumber(sample, geneCopyNumbers);
    }

    public void writeDriverCatalog(@NotNull String sample, @NotNull List<DriverCatalog> driverCatalog) {
        driverCatalogDAO.write(sample, driverCatalog);
    }

    public void writeMetrics(@NotNull String sample, @NotNull WGSMetricWithQC metrics) {
        metricDAO.writeMetrics(sample, metrics);
    }

    public void writePGX(@NotNull String sample, @NotNull List<PGXGenotype> pgxGenotype, @NotNull List<PGXCalls> pgxCalls) {
        pgxDAO.writePgx(sample, pgxGenotype, pgxCalls);
    }

    public void writeChord(@NotNull String sample, @NotNull ChordAnalysis chordAnalysis) {
        chordDAO.writeChord(sample, chordAnalysis);
    }

    public void writeClinicalEvidence(@NotNull String sample, @NotNull List<EvidenceItem> items) {
        clinicalEvidenceDAO.writeClinicalEvidence(sample, items);
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

    public void writeFullClinicalData(@NotNull Patient patient, boolean blacklisted) {
        clinicalDAO.writeFullClinicalData(patient, blacklisted);
    }

    public void writeSampleClinicalData(@NotNull String patientIdentifier, boolean blacklisted, @NotNull List<SampleData> samples) {
        clinicalDAO.writeSampleClinicalData(patientIdentifier, blacklisted, samples);
    }

    public void writeDrupEcrf(@NotNull EcrfModel model, @NotNull Set<String> sequencedPatients) {
        LOGGER.info("Writing DRUP datamodel...");
        ecrfDAO.writeDrupDatamodel(model.fields());
        LOGGER.info(" Done writing DRUP datamodel.");
        LOGGER.info("Writing raw DRUP patient data...");
        model.patients().forEach(patient -> ecrfDAO.writeDrupPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info(" Done writing raw DRUP patient data.");
    }

    public void writeCpctEcrf(@NotNull EcrfModel model, @NotNull Set<String> sequencedPatients) {
        LOGGER.info("Writing CPCT datamodel...");
        ecrfDAO.writeCpctDatamodel(model.fields());
        LOGGER.info(" Done writing CPCT datamodel.");
        LOGGER.info("Writing raw CPCT patient data...");
        model.patients().forEach(patient -> ecrfDAO.writeCpctPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info(" Done writing raw CPCT patient data.");
    }

    public void writeValidationFindings(@NotNull List<ValidationFinding> findings) {
        validationFindingsDAO.write(findings);
    }

    public void deleteAllDataForSample(@NotNull String sample) {
        LOGGER.info("Starting deleting data");

        LOGGER.info("Deleting metric data for sample: {}", sample);
        metricDAO.deleteMetricForSample(sample);

        LOGGER.info("Deleting CHORD data for sample: {}", sample);
        chordDAO.deleteChordForSample(sample);

        LOGGER.info("Deleting amber data for sample: {}", sample);
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

        LOGGER.info("Deleting structural variant annotation data for sample: {}", sample);
        structuralVariantFusionDAO.deleteAnnotationsForSample(sample);

        LOGGER.info("Deleting structural variant cluster data for sample: {}", sample);
        structuralVariantClusterDAO.deleteClusterDataForSample(sample);

        LOGGER.info("Deleting structural variants for sample: {}", sample);
        structuralVariantDAO.deleteStructuralVariantsForSample(sample);

        LOGGER.info("Deleting evidence data for sample: {}", sample);
        clinicalEvidenceDAO.deleteClinicalEvidenceForSample(sample);

        LOGGER.info("Deleting driver catalog for sample: {}", sample);
        driverCatalogDAO.deleteForSample(sample);

        LOGGER.info("Deleting pgx data for sample: {}", sample);
        pgxDAO.deletePgxForSample(sample);

        LOGGER.info("All data for sample '{}' has been deleted", sample);
    }
}

