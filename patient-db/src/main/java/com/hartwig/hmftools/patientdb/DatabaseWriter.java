package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSE;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

class DatabaseWriter {
    private static final Logger LOGGER = LogManager.getLogger(DatabaseWriter.class);

    @NotNull
    private final DSLContext context;

    DatabaseWriter(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        this.context = DSL.using(conn, SQLDialect.MYSQL);
    }

    void clearClinicalTables() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(PATIENT).execute();
        context.truncate(SAMPLE).execute();
        context.truncate(BIOPSY).execute();
        context.truncate(TREATMENT).execute();
        context.truncate(DRUG).execute();
        context.truncate(TREATMENTRESPONSE).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writeSomaticVariants(@NotNull final String sampleId,
            @NotNull final List<SomaticVariantData> somaticVariants) {
        final Record sampleRecord = context.select(SAMPLE.PATIENTID).from(SAMPLE).where(
                SAMPLE.SAMPLEID.eq(sampleId)).fetchOne();
        if (sampleRecord != null) {
            final int patientId = sampleRecord.getValue(SAMPLE.PATIENTID);
            somaticVariants.forEach(
                    somaticVariantData -> writeSomaticVariantData(patientId, sampleId, somaticVariantData));
        } else {
            LOGGER.warn(sampleId + ": was not found in table " + SAMPLE.getName());
        }
    }

    void writeClinicalData(@NotNull final Patient patient) {
        final int patientId = writePatientData(patient.patientInfo());
        patient.sequencedBiopsies().forEach(biopsy -> writeSampleData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
    }

    private int writePatientData(@NotNull final PatientData patientData) {
        final Record patientRecord = context.select(PATIENT.ID).from(PATIENT).where(
                PATIENT.CPCTID.eq(patientData.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(PATIENT.ID);
        } else {
            return context.insertInto(PATIENT, PATIENT.CPCTID, PATIENT.REGISTRATIONDATE, PATIENT.GENDER,
                    PATIENT.ETHNICITY, PATIENT.HOSPITAL, PATIENT.BIRTHYEAR, PATIENT.PRIMARYTUMORLOCATION,
                    PATIENT.DEATHDATE).values(patientData.cpctId(), Utils.toSQLDate(patientData.registrationDate()),
                    patientData.gender(), patientData.ethnicity(), patientData.hospital(), patientData.birthYear(),
                    patientData.primaryTumorLocation(), Utils.toSQLDate(patientData.deathDate())).returning(
                    PATIENT.ID).fetchOne().getValue(PATIENT.ID);
        }
    }

    private void writeSampleData(final int patientId, @NotNull final SampleData sequencedBiopsy) {
        // MIVO: ignore if primary key (sampleId) is duplicated (happens if a sample is re-sequenced).
        context.insertInto(BIOPSY, SAMPLE.SAMPLEID, SAMPLE.PATIENTID, SAMPLE.ARRIVALDATE).values(
                sequencedBiopsy.sampleId(), patientId,
                Utils.toSQLDate(sequencedBiopsy.arrivalDate())).onDuplicateKeyIgnore().execute();
    }

    private void writeBiopsyData(final int patientId, @NotNull final BiopsyClinicalData clinicalBiopsy) {
        context.insertInto(BIOPSY, BIOPSY.ID, BIOPSY.SAMPLEID, BIOPSY.PATIENTID, BIOPSY.BIOPSYLOCATION,
                BIOPSY.BIOPSYDATE).values(clinicalBiopsy.id(), clinicalBiopsy.sampleId(), patientId,
                clinicalBiopsy.location(), Utils.toSQLDate(clinicalBiopsy.date())).execute();
    }

    private void writeTreatmentData(final int patientId, @NotNull final BiopsyTreatmentData treatmentData) {
        context.insertInto(TREATMENT, TREATMENT.ID, TREATMENT.BIOPSYID, TREATMENT.PATIENTID, TREATMENT.TREATMENTGIVEN,
                TREATMENT.STARTDATE, TREATMENT.ENDDATE, TREATMENT.NAME, TREATMENT.TYPE).values(treatmentData.id(),
                treatmentData.biopsyId(), patientId, treatmentData.treatmentGiven(),
                Utils.toSQLDate(treatmentData.startDate()), Utils.toSQLDate(treatmentData.endDate()),
                treatmentData.treatmentName(), treatmentData.type()).execute();
        treatmentData.drugs().forEach(drug -> writeDrugData(patientId, treatmentData.id(), drug));
    }

    private void writeDrugData(final int patientId, final int treatmentId,
            @NotNull final BiopsyTreatmentDrugData drugData) {
        context.insertInto(DRUG, DRUG.TREATMENTID, DRUG.PATIENTID, DRUG.STARTDATE, DRUG.ENDDATE, DRUG.NAME,
                DRUG.TYPE).values(treatmentId, patientId, Utils.toSQLDate(drugData.startDate()),
                Utils.toSQLDate(drugData.endDate()), drugData.name(), drugData.type()).execute();
    }

    private void writeTreatmentResponseData(final int patientId,
            @NotNull final BiopsyTreatmentResponseData treatmentResponseData) {
        context.insertInto(TREATMENTRESPONSE, TREATMENTRESPONSE.TREATMENTID, TREATMENTRESPONSE.PATIENTID,
                TREATMENTRESPONSE.RESPONSEDATE, TREATMENTRESPONSE.RESPONSE, TREATMENTRESPONSE.MEASUREMENTDONE).values(
                treatmentResponseData.treatmentId(), patientId, Utils.toSQLDate(treatmentResponseData.date()),
                treatmentResponseData.response(), treatmentResponseData.measurementDone()).execute();
    }

    private void writeSomaticVariantData(final int patientId, @NotNull final String sampleId,
            @NotNull final SomaticVariantData somaticVariantData) {
        context.insertInto(SOMATICVARIANT, SOMATICVARIANT.SAMPLEID, SOMATICVARIANT.PATIENTID, SOMATICVARIANT.GENE,
                SOMATICVARIANT.POSITION, SOMATICVARIANT.REF, SOMATICVARIANT.ALT, SOMATICVARIANT.COSMICID,
                SOMATICVARIANT.TOTALREADCOUNT, SOMATICVARIANT.ALLELEREADCOUNT).values(sampleId, patientId,
                somaticVariantData.gene(), somaticVariantData.position(), somaticVariantData.ref(),
                somaticVariantData.alt(), somaticVariantData.cosmicID(), somaticVariantData.totalReadCount(),
                somaticVariantData.alleleReadCount()).execute();
    }
}
