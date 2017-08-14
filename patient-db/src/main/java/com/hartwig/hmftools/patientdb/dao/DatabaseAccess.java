package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.ValidationFinding;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.purity.PurityContext;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

public class DatabaseAccess {
    private static final Logger LOGGER = LogManager.getLogger(DatabaseAccess.class);

    @NotNull
    private final DSLContext context;
    @NotNull
    private final PurityDAO purityDAO;
    @NotNull
    private final CopyNumberDAO copyNumberDAO;
    @NotNull
    private final GeneCopyNumberDAO geneCopyNumberDAO;
    @NotNull
    private final ComprehensiveSomaticVariantDAO somaticVariantDAO;
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
        this.context = DSL.using(conn, SQLDialect.MYSQL);
        purityDAO = new PurityDAO(context);
        copyNumberDAO = new CopyNumberDAO(context);
        geneCopyNumberDAO = new GeneCopyNumberDAO(context);
        somaticVariantDAO = new ComprehensiveSomaticVariantDAO(context);
        structuralVariantDAO = new StructuralVariantDAO(context);
        ecrfDAO = new EcrfDAO(context);
        clinicalDAO = new ClinicalDAO(context);
        validationFindingsDAO = new ValidationFindingDAO(context);
    }

    public void writePurity(@NotNull final String sampleId, @NotNull final PurityContext context) {
        purityDAO.write(sampleId, context);
        purityDAO.write(sampleId, context.bestPerPurity());
    }

    @Nullable
    public FittedPurity readFittedPurity(@NotNull final String sampleId) {
        return purityDAO.readFittedPurity(sampleId);
    }

    @Nullable
    public FittedPurityScore readFittedPurityScore(@NotNull final String sampleId) {
        return purityDAO.readFittedPurityScore(sampleId);
    }

    @NotNull
    public List<EnrichedSomaticVariant> readComprehensiveSomaticVariants(@NotNull final String sampleId) {
        return somaticVariantDAO.read(sampleId, true);
    }

    public void writeComprehensiveSomaticVariants(@NotNull final String sampleId, @NotNull List<EnrichedSomaticVariant> variants) {
        somaticVariantDAO.write(sampleId, variants);
    }

    @NotNull
    public List<StructuralVariant> readStructuralVariants(@NotNull final String sampleId) {
        return structuralVariantDAO.read(sampleId);
    }

    public void writeStructuralVariants(@NotNull final String sampleId, @NotNull List<StructuralVariant> variants) {
        structuralVariantDAO.write(sampleId, variants);
    }

    public void writeCopynumbers(@NotNull final String sample, @NotNull List<PurpleCopyNumber> copyNumbers) {
        copyNumberDAO.writeCopyNumber(sample, copyNumbers);
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

    public void clearEcrf() {
        ecrfDAO.clear();
    }

    public void clearClinicalTables() {
        validationFindingsDAO.clear();
        clinicalDAO.clear();
    }

    public void clearSomaticTables() {
        context.truncate(SOMATICVARIANT).execute();
    }

    public void writeClinicalData(@NotNull final Patient patient) {
        clinicalDAO.writeClinicalData(patient);
    }

    public void writeSomaticVariants(@NotNull final String sampleId, @NotNull final List<SomaticVariantData> somaticVariants) {
        context.batch(somaticVariants.stream()
                .map(somaticVariant -> context.insertInto(SOMATICVARIANT, SOMATICVARIANT.SAMPLEID, SOMATICVARIANT.GENE,
                        SOMATICVARIANT.POSITION, SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.COSMICID,
                        SOMATICVARIANT.TOTALREADCOUNT, SOMATICVARIANT.ALLELEREADCOUNT)
                        .values(sampleId, somaticVariant.gene(), somaticVariant.position(), somaticVariant.ref(), somaticVariant.alt(),
                                somaticVariant.cosmicID(), somaticVariant.totalReadCount(), somaticVariant.alleleReadCount()))
                .collect(Collectors.toList())).execute();
    }

    public void writeEcrf(@NotNull final CpctEcrfModel model, @NotNull final Set<String> sequencedPatients) {
        LOGGER.info("writing datamodel...");
        ecrfDAO.writeDatamodel(model.fields());
        LOGGER.info("done writing datamodel.");
        LOGGER.info("writing patients...");
        model.patients().forEach(patient -> ecrfDAO.writePatient(patient, sequencedPatients.contains(patient.patientId())));
        LOGGER.info("done writing patients.");
    }

    public void writeValidationFindings(@NotNull final List<ValidationFinding> findings) {
        validationFindingsDAO.write(findings);
    }
}
