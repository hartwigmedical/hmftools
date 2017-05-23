package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.CLINICALBIOPSIES;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUGS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.LIMSBIOPSIES;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENTS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANTS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSES;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTS;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;
import java.util.List;

import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientInfo;
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
        context.truncate(PATIENTS).execute();
        context.truncate(LIMSBIOPSIES).execute();
        context.truncate(CLINICALBIOPSIES).execute();
        context.truncate(TREATMENTS).execute();
        context.truncate(DRUGS).execute();
        context.truncate(TREATMENTRESPONSES).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writeSomaticVariants(@NotNull final String sampleId,
            @NotNull final List<SomaticVariantData> somaticVariants) {
        final Record limsRecord = context.select(LIMSBIOPSIES.PATIENTID).from(LIMSBIOPSIES).where(
                LIMSBIOPSIES.SAMPLEID.eq(sampleId)).fetchOne();
        if (limsRecord != null) {
            final int patientId = limsRecord.getValue(LIMSBIOPSIES.PATIENTID);
            somaticVariants.forEach(
                    somaticVariantData -> writeSomaticVariantData(patientId, sampleId, somaticVariantData));
        } else {
            LOGGER.warn(sampleId + ": was not found in table " + LIMSBIOPSIES.getName());
        }
    }

    void writeClinicalData(@NotNull final Patient patient) {
        final int patientId = writePatientInfo(patient.patientInfo());
        patient.sequencedBiopsies().forEach(biopsy -> writeBiopsyLimsData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeClinicalBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
    }

    private int writePatientInfo(@NotNull final PatientInfo patientInfo) {
        final Record patientRecord = context.select(PATIENTS.ID).from(PATIENTS).where(
                PATIENTS.CPCTID.eq(patientInfo.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(PATIENTS.ID);
        } else {
            return context.insertInto(PATIENTS, PATIENTS.CPCTID, PATIENTS.REGISTRATIONDATE, PATIENTS.GENDER,
                    PATIENTS.ETHNICITY, PATIENTS.HOSPITAL, PATIENTS.BIRTHYEAR, PATIENTS.TUMORLOCATION,
                    PATIENTS.DEATHDATE).values(patientInfo.cpctId(), Utils.toSQLDate(patientInfo.registrationDate()),
                    patientInfo.gender(), patientInfo.ethnicity(), patientInfo.hospital(), patientInfo.birthYear(),
                    patientInfo.primaryTumorLocation(), Utils.toSQLDate(patientInfo.deathDate())).returning(
                    PATIENTS.ID).fetchOne().getValue(PATIENTS.ID);
        }
    }

    private void writeBiopsyLimsData(final int patientId, @NotNull final BiopsyLimsData sequencedBiopsy) {
        // MIVO: ignore if primary key (sampleId) is duplicated (happens if a sample is re-sequenced).
        context.insertInto(LIMSBIOPSIES, LIMSBIOPSIES.SAMPLEID, LIMSBIOPSIES.ARRIVALDATE,
                LIMSBIOPSIES.PATIENTID).values(sequencedBiopsy.sampleId(),
                Utils.toSQLDate(sequencedBiopsy.arrivalDate()), patientId).onDuplicateKeyIgnore().execute();
    }

    private void writeClinicalBiopsyData(final int patientId, @NotNull final BiopsyClinicalData clinicalBiopsy) {
        context.insertInto(CLINICALBIOPSIES, CLINICALBIOPSIES.ID, CLINICALBIOPSIES.LOCATION, CLINICALBIOPSIES.DATE,
                CLINICALBIOPSIES.SAMPLEID, CLINICALBIOPSIES.PATIENTID).values(clinicalBiopsy.id(),
                clinicalBiopsy.location(), Utils.toSQLDate(clinicalBiopsy.date()), clinicalBiopsy.sampleId(),
                patientId).execute();
    }

    private void writeTreatmentData(final int patientId, @NotNull final BiopsyTreatmentData treatmentData) {
        context.insertInto(TREATMENTS, TREATMENTS.ID, TREATMENTS.TREATMENTGIVEN, TREATMENTS.STARTDATE,
                TREATMENTS.ENDDATE, TREATMENTS.NAME, TREATMENTS.TYPE, TREATMENTS.CLINICALBIOPSYID,
                TREATMENTS.PATIENTID).values(treatmentData.id(), treatmentData.treatmentGiven(),
                Utils.toSQLDate(treatmentData.startDate()), Utils.toSQLDate(treatmentData.endDate()),
                treatmentData.treatmentName(), treatmentData.type(), treatmentData.biopsyId(), patientId).execute();
        treatmentData.drugs().forEach(drug -> writeDrugData(patientId, treatmentData.id(), drug));
    }

    private void writeDrugData(final int patientId, final int treatmentId,
            @NotNull final BiopsyTreatmentDrugData drugData) {
        context.insertInto(DRUGS, DRUGS.STARTDATE, DRUGS.ENDDATE, DRUGS.NAME, DRUGS.TYPE, DRUGS.TREATMENTID,
                DRUGS.PATIENTID).values(Utils.toSQLDate(drugData.startDate()), Utils.toSQLDate(drugData.endDate()),
                drugData.name(), drugData.type(), treatmentId, patientId).execute();
    }

    private void writeTreatmentResponseData(final int patientId,
            @NotNull final BiopsyTreatmentResponseData treatmentResponseData) {
        context.insertInto(TREATMENTRESPONSES, TREATMENTRESPONSES.DATE, TREATMENTRESPONSES.RESPONSE,
                TREATMENTRESPONSES.MEASUREMENTDONE, TREATMENTRESPONSES.TREATMENTID,
                TREATMENTRESPONSES.PATIENTID).values(Utils.toSQLDate(treatmentResponseData.date()),
                treatmentResponseData.response(), treatmentResponseData.measurementDone(),
                treatmentResponseData.treatmentId(), patientId).execute();
    }

    private void writeSomaticVariantData(final int patientId, @NotNull final String sampleId,
            @NotNull final SomaticVariantData somaticVariantData) {
        context.insertInto(SOMATICVARIANTS, SOMATICVARIANTS.GENE, SOMATICVARIANTS.POSITION, SOMATICVARIANTS.REF,
                SOMATICVARIANTS.ALT, SOMATICVARIANTS.COSMICID, SOMATICVARIANTS.TOTALREADCOUNT,
                SOMATICVARIANTS.ALLELEREADCOUNT, SOMATICVARIANTS.PATIENTID, SOMATICVARIANTS.SAMPLEID).values(
                somaticVariantData.gene(), somaticVariantData.position(), somaticVariantData.ref(),
                somaticVariantData.alt(), somaticVariantData.cosmicID(), somaticVariantData.totalReadCount(),
                somaticVariantData.alleleReadCount(), patientId, sampleId).execute();
    }
}
