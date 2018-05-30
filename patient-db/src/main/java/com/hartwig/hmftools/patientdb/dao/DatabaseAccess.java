package com.hartwig.hmftools.patientdb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.common.ecrf.EcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.gene.CanonicalTranscript;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.metrics.WGSMetrics;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PotentialActionableCNV;
import com.hartwig.hmftools.patientdb.data.PotentialActionableFusion;
import com.hartwig.hmftools.patientdb.data.PotentialActionableVariant;
import com.hartwig.hmftools.patientdb.data.SampleData;

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

public class DatabaseAccess {
    private static final Logger LOGGER = LogManager.getLogger(DatabaseAccess.class);
    private static final String DEV_CATALOG = "hmfpatients";

    @NotNull
    private final DSLContext context;
    @NotNull
    private final PurityDAO purityDAO;
    @NotNull
    private final CopyNumberDAO copyNumberDAO;
    @NotNull
    private final GeneCopyNumberDAO geneCopyNumberDAO;
    @NotNull
    private final SomaticVariantDAO somaticVariantDAO;
    @NotNull
    private final StructuralVariantDAO structuralVariantDAO;
    @NotNull
    private final EcrfDAO ecrfDAO;
    @NotNull
    private final ClinicalDAO clinicalDAO;
    @NotNull
    private final ValidationFindingDAO validationFindingsDAO;
    @NotNull
    private final CanonicalTranscriptDAO canonicalTranscriptDAO;
    @NotNull
    private final MetricDAO metricDAO;
    @NotNull
    private final PotentiallyActionableItemsDAO potentiallyActionableItemsDAO;

    public DatabaseAccess(@NotNull final String userName, @NotNull final String password, @NotNull final String url) throws SQLException {
        // MIVO: disable annoying jooq self-ad message
        System.setProperty("org.jooq.no-logo", "true");
        final Connection conn = DriverManager.getConnection(url, userName, password);
        final String catalog = conn.getCatalog();
        LOGGER.info("Connecting to database {}", catalog);
        this.context = DSL.using(conn, SQLDialect.MYSQL, settings(catalog));

        purityDAO = new PurityDAO(context);
        copyNumberDAO = new CopyNumberDAO(context);
        geneCopyNumberDAO = new GeneCopyNumberDAO(context);
        somaticVariantDAO = new SomaticVariantDAO(context);
        structuralVariantDAO = new StructuralVariantDAO(context);
        ecrfDAO = new EcrfDAO(context);
        clinicalDAO = new ClinicalDAO(context);
        validationFindingsDAO = new ValidationFindingDAO(context);
        canonicalTranscriptDAO = new CanonicalTranscriptDAO(context);
        metricDAO = new MetricDAO(context);
        potentiallyActionableItemsDAO = new PotentiallyActionableItemsDAO(context);
    }

    @NotNull
    public DSLContext context() {
        return context;
    }

    @Nullable
    private static Settings settings(final String catalog) {
        return !catalog.equals(DEV_CATALOG)
                ? new Settings().withRenderMapping(new RenderMapping().withSchemata(new MappedSchema().withInput(DEV_CATALOG)
                .withOutput(catalog)))
                : null;
    }

    public void writeCanonicalTranscripts(@NotNull final List<CanonicalTranscript> transcripts) {
        canonicalTranscriptDAO.write(transcripts);
    }

    public void writePurity(@NotNull final String sampleId, @NotNull final PurityContext context, @NotNull final PurpleQC checks) {
        purityDAO.write(sampleId, context, checks);
    }

    public void writeBestFitPerPurity(@NotNull final String sampleId, @NotNull final List<FittedPurity> bestFitPerPurity) {
        purityDAO.write(sampleId, bestFitPerPurity);
    }

    @Nullable
    public PurityContext readPurityContext(@NotNull final String sampleId) {
        return purityDAO.readPurityContext(sampleId);
    }

    public void writeSomaticVariants(@NotNull final String sampleId, @NotNull List<EnrichedSomaticVariant> variants) {
        somaticVariantDAO.write(sampleId, variants);
    }

    @NotNull
    public Stream<PotentialActionableVariant> potentiallyActionableVariants() {
        return potentiallyActionableItemsDAO.potentiallyActionableVariants();
    }

    @NotNull
    public Stream<PotentialActionableCNV> potentiallyActionableCNVs() {
        return potentiallyActionableItemsDAO.potentiallyActionableCNVs();
    }

    @NotNull
    public Stream<PotentialActionableFusion> potentiallyActionableFusions() {
        return potentiallyActionableItemsDAO.potentiallyActionableFusions();
    }

    public void writeStructuralVariants(@NotNull final String sampleId, @NotNull final List<EnrichedStructuralVariant> variants) {
        structuralVariantDAO.write(sampleId, variants);
    }

    public void writeCopynumbers(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeCopyNumber(sample, copyNumbers);
    }

    @NotNull
    public List<StructuralVariantData> readStructuralVariantData(@NotNull final String sample) {
        return structuralVariantDAO.read(sample);
    }

    @NotNull
    public List<EnrichedStructuralVariant> readStructuralVariants(@NotNull final String sample) {
        return structuralVariantDAO.readEnrichedData(sample);
    }

    @NotNull
    public List<String> structuralVariantSampleList(@NotNull final String sampleSearch) {
        return structuralVariantDAO.getSamplesList(sampleSearch);
    }

    public void writeGermlineCopynumbers(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeGermlineCopyNumber(sample, copyNumbers);
    }

    public void writeCopynumberRegions(@NotNull final String sample, @NotNull List<FittedRegion> regions) {
        copyNumberDAO.writeCopyNumberRegions(sample, regions);
    }

    @NotNull
    public List<FittedRegion> readCopyNumberRegions(@NotNull final String sample) {
        return copyNumberDAO.readCopyNumberRegions(sample);
    }

    public void writeGeneCopynumberRegions(@NotNull final String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        geneCopyNumberDAO.writeCopyNumber(sample, geneCopyNumbers);
    }

    @NotNull
    public List<PurpleCopyNumber> readCopynumbers(@NotNull final String sample) {
        return copyNumberDAO.read(sample);
    }

    public void writeMetrics(@NotNull String sample, @NotNull WGSMetrics metrics) {
        metricDAO.writeMetrics(sample, metrics);
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

    public void writeFullClinicalData(@NotNull Patient patient) {
        clinicalDAO.writeFullClinicalData(patient);
    }

    public void writeSampleClinicalData(@NotNull String patientIdentifier, @NotNull List<SampleData> samples) {
        clinicalDAO.writeSampleClinicalData(patientIdentifier, samples);
    }

    public void writeDrupEcrf(@NotNull final EcrfModel model, @NotNull final Set<String> sequencedPatients) {
        LOGGER.info("Writing DRUP datamodel...");
        ecrfDAO.writeDrupDatamodel(model.fields());
        LOGGER.info("Done writing DRUP datamodel.");
        LOGGER.info("Writing DRUP patients...");
        model.patients().forEach(patient -> ecrfDAO.writeDrupPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("Done writing DRUP patients.");
    }

    public void writeCpctEcrf(@NotNull final EcrfModel model, @NotNull final Set<String> sequencedPatients) {
        LOGGER.info("Writing CPCT datamodel...");
        ecrfDAO.writeCpctDatamodel(model.fields());
        LOGGER.info("Done writing CPCT datamodel.");
        LOGGER.info("Writing CPCT patients...");
        model.patients().forEach(patient -> ecrfDAO.writeCpctPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("Done writing CPCT patients.");
    }

    public void writeValidationFindings(@NotNull final List<ValidationFinding> findings) {
        validationFindingsDAO.write(findings);
    }

    public void deleteAllDataForSample(@NotNull String sample) {
        LOGGER.info("Starting deleting data");

        LOGGER.info("Deleting metric data for sample: " + sample);
        metricDAO.deleteMetricForSample(sample);

        LOGGER.info("Deleting purity data for sample: " + sample);
        purityDAO.deletePurityForSample(sample);

        LOGGER.info("Deleting copy number data for sample: " + sample);
        copyNumberDAO.deleteCopyNumberForSample(sample);

        LOGGER.info("Deleting gene copy number data for sample: " + sample);
        geneCopyNumberDAO.deleteGeneCopyNumberForSample(sample);

        LOGGER.info("Deleting somatic variant data for sample: " + sample);
        somaticVariantDAO.deleteSomaticVariantForSample(sample);

        LOGGER.info("Deleting structural variant data for sample: " + sample);
        structuralVariantDAO.deleteStructuralVariantsForSample(sample);

        LOGGER.info("TODO: remove set from datastore");
        LOGGER.info("TODO: run run_clinical_patient_db_prod");
        LOGGER.info("Clinical data for sample: " + sample + " is removed from database");

        LOGGER.info("All data for sample: " + sample + " is deleted");

    }
}

