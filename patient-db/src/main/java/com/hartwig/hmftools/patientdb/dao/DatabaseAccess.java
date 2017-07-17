package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALLOGS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.FORMSMETADATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SEQUENCEDPATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSE;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.gene.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.copynumber.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.purity.FittedPurity;
import com.hartwig.hmftools.common.purple.purity.FittedPurityScore;
import com.hartwig.hmftools.common.purple.region.FittedRegion;
import com.hartwig.hmftools.common.variant.EnrichedSomaticVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.jooq.DSLContext;
import org.jooq.Record;
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
    }

    public void writePurity(@NotNull final String sampleId, @NotNull FittedPurityScore score, @NotNull List<FittedPurity> purities) {
        purityDAO.write(sampleId, score);
        purityDAO.write(sampleId, purities);
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

    public void clearClinicalTables() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(CLINICALLOGS).execute();
        context.truncate(PATIENT).execute();
        context.truncate(SAMPLE).execute();
        context.truncate(BIOPSY).execute();
        context.truncate(TREATMENT).execute();
        context.truncate(DRUG).execute();
        context.truncate(TREATMENTRESPONSE).execute();
        context.truncate(FORMSMETADATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    public void clearEcrf() {
        ecrfDAO.clear();
    }

    public void writeSomaticVariants(@NotNull final String sampleId, @NotNull final List<SomaticVariantData> somaticVariants) {
        final Record sampleRecord = context.select(SAMPLE.PATIENTID).from(SAMPLE).where(SAMPLE.SAMPLEID.eq(sampleId)).fetchOne();
        if (sampleRecord != null) {
            final int patientId = sampleRecord.getValue(SAMPLE.PATIENTID);
            somaticVariants.forEach(somaticVariantData -> writeSomaticVariantData(patientId, sampleId, somaticVariantData));
        } else {
            LOGGER.warn(sampleId + ": was not found in table " + SAMPLE.getName());
        }
    }

    public void writeClinicalData(@NotNull final Patient patient) {
        final int patientId = writeSequencedPatient(patient.patientInfo());
        writePatientData(patientId, patient.patientInfo());
        patient.sequencedBiopsies().forEach(biopsy -> writeSampleData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
    }

    private int writeSequencedPatient(@NotNull final PatientData patient) {
        final Record patientRecord =
                context.select(SEQUENCEDPATIENT.ID).from(SEQUENCEDPATIENT).where(SEQUENCEDPATIENT.CPCTID.eq(patient.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(SEQUENCEDPATIENT.ID);
        } else {
            return context.insertInto(SEQUENCEDPATIENT, SEQUENCEDPATIENT.CPCTID)
                    .values(patient.cpctId())
                    .returning(SEQUENCEDPATIENT.ID)
                    .fetchOne()
                    .getValue(SEQUENCEDPATIENT.ID);

        }
    }

    private void writePatientData(final int patientId, @NotNull final PatientData patient) {
        context.insertInto(PATIENT, PATIENT.ID, PATIENT.REGISTRATIONDATE, PATIENT.GENDER, PATIENT.ETHNICITY, PATIENT.HOSPITAL,
                PATIENT.BIRTHYEAR, PATIENT.PRIMARYTUMORLOCATION, PATIENT.DEATHDATE)
                .values(patientId, Utils.toSQLDate(patient.registrationDate()), patient.gender(), patient.ethnicity(), patient.hospital(),
                        patient.birthYear(), patient.primaryTumorLocation(), Utils.toSQLDate(patient.deathDate()))
                .execute();
        writeFormStatus(patientId, PATIENT.getName(), "demography", patient.demographyStatus(), patient.demographyLocked());
        writeFormStatus(patientId, PATIENT.getName(), "primaryTumor", patient.primaryTumorStatus(), patient.primaryTumorLocked());
        writeFormStatus(patientId, PATIENT.getName(), "eligibility", patient.eligibilityStatus(), patient.eligibilityLocked());
        writeFormStatus(patientId, PATIENT.getName(), "selectionCriteria", patient.selectionCriteriaStatus(),
                patient.selectionCriteriaLocked());
        writeFormStatus(patientId, PATIENT.getName(), "death", patient.deathStatus(), patient.deathLocked());
    }

    private void writeSampleData(final int patientId, @NotNull final SampleData sample) {
        // MIVO: ignore if primary key (sampleId) is duplicated (happens if a sample is re-sequenced).
        context.insertInto(SAMPLE, SAMPLE.SAMPLEID, SAMPLE.PATIENTID, SAMPLE.ARRIVALDATE)
                .values(sample.sampleId(), patientId, Utils.toSQLDate(sample.arrivalDate()))
                .onDuplicateKeyIgnore()
                .execute();
    }

    private void writeBiopsyData(final int patientId, @NotNull final BiopsyData biopsy) {
        context.insertInto(BIOPSY, BIOPSY.ID, BIOPSY.SAMPLEID, BIOPSY.PATIENTID, BIOPSY.BIOPSYLOCATION, BIOPSY.BIOPSYDATE)
                .values(biopsy.id(), biopsy.sampleId(), patientId, biopsy.location(), Utils.toSQLDate(biopsy.date()))
                .execute();
        writeFormStatus(biopsy.id(), BIOPSY.getName(), "biopsy", biopsy.formStatus(), biopsy.formLocked());
    }

    private void writeTreatmentData(final int patientId, @NotNull final BiopsyTreatmentData treatment) {
        context.insertInto(TREATMENT, TREATMENT.ID, TREATMENT.BIOPSYID, TREATMENT.PATIENTID, TREATMENT.TREATMENTGIVEN, TREATMENT.STARTDATE,
                TREATMENT.ENDDATE, TREATMENT.NAME, TREATMENT.TYPE)
                .values(treatment.id(), treatment.biopsyId(), patientId, treatment.treatmentGiven(), Utils.toSQLDate(treatment.startDate()),
                        Utils.toSQLDate(treatment.endDate()), treatment.treatmentName(), treatment.type())
                .execute();
        writeFormStatus(treatment.id(), TREATMENT.getName(), "treatment", treatment.formStatus(), treatment.formLocked());
        treatment.drugs().forEach(drug -> writeDrugData(patientId, treatment.id(), drug, treatment.formStatus(), treatment.formLocked()));
    }

    private void writeDrugData(final int patientId, final int treatmentId, @NotNull final BiopsyTreatmentDrugData drug,
            @NotNull final String formStatus, @NotNull final String formLocked) {
        final int id = context.insertInto(DRUG, DRUG.TREATMENTID, DRUG.PATIENTID, DRUG.STARTDATE, DRUG.ENDDATE, DRUG.NAME, DRUG.TYPE)
                .values(treatmentId, patientId, Utils.toSQLDate(drug.startDate()), Utils.toSQLDate(drug.endDate()), drug.name(),
                        drug.type())
                .returning(DRUG.ID)
                .fetchOne()
                .getValue(DRUG.ID);
        writeFormStatus(id, DRUG.getName(), "treatment", formStatus, formLocked);
    }

    private void writeTreatmentResponseData(final int patientId, @NotNull final BiopsyTreatmentResponseData treatmentResponse) {
        final int id = context.insertInto(TREATMENTRESPONSE, TREATMENTRESPONSE.TREATMENTID, TREATMENTRESPONSE.PATIENTID,
                TREATMENTRESPONSE.RESPONSEDATE, TREATMENTRESPONSE.RESPONSE, TREATMENTRESPONSE.MEASUREMENTDONE)
                .values(treatmentResponse.treatmentId(), patientId, Utils.toSQLDate(treatmentResponse.date()), treatmentResponse.response(),
                        treatmentResponse.measurementDone())
                .returning(TREATMENTRESPONSE.ID)
                .fetchOne()
                .getValue(TREATMENTRESPONSE.ID);
        writeFormStatus(id, TREATMENTRESPONSE.getName(), "treatmentResponse", treatmentResponse.formStatus(),
                treatmentResponse.formLocked());
    }

    private void writeSomaticVariantData(final int patientId, @NotNull final String sampleId,
            @NotNull final SomaticVariantData somaticVariant) {
        context.insertInto(SOMATICVARIANT, SOMATICVARIANT.SAMPLEID, SOMATICVARIANT.PATIENTID, SOMATICVARIANT.GENE, SOMATICVARIANT.POSITION,
                SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.COSMICID, SOMATICVARIANT.TOTALREADCOUNT,
                SOMATICVARIANT.ALLELEREADCOUNT)
                .values(sampleId, patientId, somaticVariant.gene(), somaticVariant.position(), somaticVariant.ref(), somaticVariant.alt(),
                        somaticVariant.cosmicID(), somaticVariant.totalReadCount(), somaticVariant.alleleReadCount())
                .execute();
    }

    private void writeFormStatus(final int id, @NotNull final String tableName, @NotNull final String formName,
            @NotNull final String formStatus, @NotNull final String formLocked) {
        context.insertInto(FORMSMETADATA, FORMSMETADATA.ID, FORMSMETADATA.TABLENAME, FORMSMETADATA.FORM, FORMSMETADATA.STATUS,
                FORMSMETADATA.LOCKED).values(id, tableName, formName, formStatus, formLocked).execute();
    }

    public void writeEcrfPatients(@NotNull final Iterable<EcrfPatient> patients, @NotNull final Set<String> sequencedPatients) {
        patients.forEach(patient -> ecrfDAO.writePatient(patient, sequencedPatients.contains(patient.patientId())));
    }
}
