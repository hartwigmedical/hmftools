package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSYLOCATIONS;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENTINFO;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.RADIOTHERAPYDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SYSTEMICTHERAPYDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TUMORDATA;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.SQLDialect;
import org.jooq.impl.DSL;

class DatabaseWriter {
    @NotNull
    private final DSLContext context;

    DatabaseWriter(@NotNull final String userName, @NotNull final String password, @NotNull final String url)
            throws SQLException {
        final Connection conn = DriverManager.getConnection(url, userName, password);
        this.context = DSL.using(conn, SQLDialect.MYSQL);
    }

    void clearTables() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(PATIENTINFO).execute();
        context.truncate(TUMORDATA).execute();
        context.truncate(BIOPSYLOCATIONS).execute();
        context.truncate(TREATMENTDATA).execute();
        context.truncate(RADIOTHERAPYDATA).execute();
        context.truncate(SYSTEMICTHERAPYDATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writePatient(@NotNull final Patient patient) {
        final int patientId = writePatientInfo(patient.patientInfo());
        patient.tumorData().ifPresent(tumorData -> writeTumorData(patientId, tumorData));
        patient.systemicTherapies().ifPresent(systemicTherapies -> systemicTherapies.forEach(
                systemicTherapyData -> writeSystemicTherapyData(patientId, systemicTherapyData)));
        patient.radioTherapies().ifPresent(radioTherapies -> radioTherapies.forEach(
                radioTherapyData -> writeRadioTherapyData(patientId, radioTherapyData)));
        patient.treatmentData().ifPresent(treatmentData -> writeTreatmentData(patientId, treatmentData));
    }

    private int writePatientInfo(@NotNull final PatientInfo patientInfo) {
        return context.insertInto(PATIENTINFO, PATIENTINFO.CPCTID, PATIENTINFO.SEX, PATIENTINFO.BIRTHYEAR,
                PATIENTINFO.HOSPITAL, PATIENTINFO.ETHNICITY).values(patientInfo.cpctId(), patientInfo.sex(),
                patientInfo.birthYear(), patientInfo.hospital(), patientInfo.ethnicity()).returning(
                PATIENTINFO.ID).fetchOne().getValue(PATIENTINFO.ID);
    }

    private void writeTumorData(final int patientId, @NotNull final TumorData tumorData) {
        final int tumorDataId = context.insertInto(TUMORDATA, TUMORDATA.LOCATION, TUMORDATA.ENTRYSTAGE,
                TUMORDATA.PATIENTID).values(tumorData.location(), tumorData.entryStage(), patientId).returning(
                TUMORDATA.ID).fetchOne().getValue(TUMORDATA.ID);
        for (final String biopsyLocation : tumorData.biopsyLocations()) {
            context.insertInto(BIOPSYLOCATIONS, BIOPSYLOCATIONS.LOCATION, BIOPSYLOCATIONS.TUMORID).values(
                    biopsyLocation, tumorDataId).execute();
        }
    }

    private void writeTreatmentData(final int patientId, @NotNull final TreatmentData treatmentData) {
        context.insertInto(TREATMENTDATA, TREATMENTDATA.NAME, TREATMENTDATA.STARTDATE, TREATMENTDATA.ENDDATE,
                TREATMENTDATA.EARLYRESPONSE, TREATMENTDATA.PATIENTID).values(treatmentData.treatmentName(),
                Utils.toSQLDate(treatmentData.startDate()), Utils.toSQLDate(treatmentData.endDate()),
                treatmentData.earlyResponse(), patientId).execute();
    }

    private void writeRadioTherapyData(final int patientId, @NotNull final RadioTherapyData radioTherapyData) {
        context.insertInto(RADIOTHERAPYDATA, RADIOTHERAPYDATA.SITE, RADIOTHERAPYDATA.ENDDATE,
                RADIOTHERAPYDATA.PATIENTID).values(radioTherapyData.site(),
                Utils.toSQLDate(radioTherapyData.endDate()), patientId).execute();
    }

    private void writeSystemicTherapyData(final int patientId,
            @NotNull final SystemicTherapyData systemicTherapyData) {
        context.insertInto(SYSTEMICTHERAPYDATA, SYSTEMICTHERAPYDATA.STARTDATE, SYSTEMICTHERAPYDATA.ENDDATE,
                SYSTEMICTHERAPYDATA.TYPE, SYSTEMICTHERAPYDATA.TREATMENT, SYSTEMICTHERAPYDATA.BESTRESPONSE,
                SYSTEMICTHERAPYDATA.PATIENTID).values(Utils.toSQLDate(systemicTherapyData.startDate()),
                Utils.toSQLDate(systemicTherapyData.endDate()), systemicTherapyData.type(),
                systemicTherapyData.treatment(), systemicTherapyData.bestResponse(), patientId).execute();
    }
}
