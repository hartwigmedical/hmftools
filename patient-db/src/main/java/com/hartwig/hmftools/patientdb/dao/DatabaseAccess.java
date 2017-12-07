package com.hartwig.hmftools.patientdb.dao;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.qc.PurpleQC;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.EnrichedStructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientdb.data.Patient;

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
    }

    @NotNull
    public DSLContext getContext() {
        return context;
    }

    @Nullable
    private Settings settings(final String catalog) {
        return !catalog.equals(DEV_CATALOG) ? new Settings().withRenderMapping(
                new RenderMapping().withSchemata(new MappedSchema().withInput("hmfpatients").withOutput(catalog))) : null;
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

    public void writeStructuralVariants(@NotNull final String sampleId, @NotNull final List<EnrichedStructuralVariant> variants) {
        structuralVariantDAO.write(sampleId, variants);
    }

    @NotNull
    public List<StructuralVariant> readStructuralVariants(@NotNull final String sample) {
        return structuralVariantDAO.read(sample);
    }

    public void writeCopynumbers(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeCopyNumber(sample, copyNumbers);
    }

    public void writeGermlineCopynumbers(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeGermlineCopyNumber(sample, copyNumbers);
    }

    public void writeCopynumberRegions(@NotNull final String sample, @NotNull List<FittedRegion> regions) {
        copyNumberDAO.writeCopyNumberRegions(sample, regions);
    }

    public void writeGeneCopynumberRegions(@NotNull final String sample, @NotNull List<GeneCopyNumber> geneCopyNumbers) {
        geneCopyNumberDAO.writeCopyNumber(sample, geneCopyNumbers);
    }

    @NotNull
    public List<PurpleCopyNumber> readCopynumbers(@NotNull final String sample) {
        return copyNumberDAO.read(sample);
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

    public void writeClinicalData(@NotNull final Patient patient) {
        clinicalDAO.writeClinicalData(patient);
    }

    public void writeDrupEcrf(@NotNull final CpctEcrfModel model, @NotNull final Set<String> sequencedPatients) {
        LOGGER.info("writing DRUP datamodel...");
        ecrfDAO.writeDrupDatamodel(model.fields());
        LOGGER.info("done writing DRUP datamodel.");
        LOGGER.info("writing DRUP patients...");
        model.patients().forEach(patient -> ecrfDAO.writeDrupPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("done writing DRUP patients.");
    }

    public void writeCpctEcrf(@NotNull final CpctEcrfModel model, @NotNull final Set<String> sequencedPatients) {
        LOGGER.info("writing CPCT datamodel...");
        ecrfDAO.writeCpctDatamodel(model.fields());
        LOGGER.info("done writing CPCT datamodel.");
        LOGGER.info("writing CPCT patients...");
        model.patients().forEach(patient -> ecrfDAO.writeCpctPatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("done writing CPCT patients.");
    }

    public void writeValidationFindings(@NotNull final List<ValidationFinding> findings) {
        validationFindingsDAO.write(findings);
    }
}
