package com.hartwig.hmftools.patientdb;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSYCLINICALDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSYLIMSDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUGDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENTINFO;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SOMATICVARIANTDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTDATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSEDATA;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.SQLException;

import com.hartwig.hmftools.patientdb.data.BiopsyClinicalData;
import com.hartwig.hmftools.patientdb.data.BiopsyLimsData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientInfo;
import com.hartwig.hmftools.patientdb.data.SomaticVariantData;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;
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
        context.truncate(BIOPSYLIMSDATA).execute();
        context.truncate(BIOPSYCLINICALDATA).execute();
        context.truncate(TREATMENTDATA).execute();
        context.truncate(DRUGDATA).execute();
        context.truncate(TREATMENTRESPONSEDATA).execute();
        context.truncate(SOMATICVARIANTDATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writePatient(@NotNull final Patient patient) {
        final int patientId = writePatientInfo(patient.patientInfo());
        patient.sequencedBiopsies().forEach(biopsy -> writeBiopsyLimsData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeClinicalBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
        patient.somaticVariants().forEach(
                somaticVariantData -> writeSomaticVariantData(patientId, somaticVariantData));
    }

    private int writePatientInfo(@NotNull final PatientInfo patientInfo) {
        final Record patientRecord = context.select(PATIENTINFO.ID).from(PATIENTINFO).where(
                PATIENTINFO.CPCTID.eq(patientInfo.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(PATIENTINFO.ID);
        } else {
            return context.insertInto(PATIENTINFO, PATIENTINFO.CPCTID, PATIENTINFO.REGISTRATIONDATE,
                    PATIENTINFO.GENDER, PATIENTINFO.ETHNICITY, PATIENTINFO.HOSPITAL, PATIENTINFO.BIRTHYEAR,
                    PATIENTINFO.TUMORLOCATION, PATIENTINFO.DEATHDATE).values(patientInfo.cpctId(),
                    Utils.toSQLDate(patientInfo.registrationDate()), patientInfo.gender(), patientInfo.ethnicity(),
                    patientInfo.hospital(), patientInfo.birthYear(), patientInfo.tumorLocation(),
                    Utils.toSQLDate(patientInfo.deathDate())).returning(PATIENTINFO.ID).fetchOne().getValue(
                    PATIENTINFO.ID);
        }
    }

    private void writeBiopsyLimsData(final int patientId, @NotNull final BiopsyLimsData sequencedBiopsy) {
        // MIVO: ignore if primary key (sampleId) is duplicated (happens if a sample is re-sequenced).
        context.insertInto(BIOPSYLIMSDATA, BIOPSYLIMSDATA.SAMPLEID, BIOPSYLIMSDATA.ARRIVALDATE,
                BIOPSYLIMSDATA.PATIENTID).values(sequencedBiopsy.sampleId(),
                Utils.toSQLDate(sequencedBiopsy.arrivalDate()), patientId).onDuplicateKeyIgnore().execute();
    }

    private void writeClinicalBiopsyData(final int patientId, @NotNull final BiopsyClinicalData clinicalBiopsy) {
        context.insertInto(BIOPSYCLINICALDATA, BIOPSYCLINICALDATA.ID, BIOPSYCLINICALDATA.LOCATION,
                BIOPSYCLINICALDATA.DATE, BIOPSYCLINICALDATA.SAMPLEID, BIOPSYCLINICALDATA.PATIENTID).values(
                clinicalBiopsy.id(), clinicalBiopsy.location(), Utils.toSQLDate(clinicalBiopsy.date()),
                clinicalBiopsy.sampleId(), patientId).execute();
    }

    private void writeTreatmentData(final int patientId, @NotNull final BiopsyTreatmentData treatmentData) {
        context.insertInto(TREATMENTDATA, TREATMENTDATA.ID, TREATMENTDATA.TREATMENTGIVEN, TREATMENTDATA.STARTDATE,
                TREATMENTDATA.ENDDATE, TREATMENTDATA.NAME, TREATMENTDATA.TYPE, TREATMENTDATA.CLINICALBIOPSYID,
                TREATMENTDATA.PATIENTID).values(treatmentData.id(), treatmentData.treatmentGiven(),
                Utils.toSQLDate(treatmentData.startDate()), Utils.toSQLDate(treatmentData.endDate()),
                treatmentData.treatmentName(), treatmentData.type(), treatmentData.biopsyId(), patientId).execute();
        treatmentData.drugs().forEach(drug -> writeDrugData(patientId, treatmentData.id(), drug));
    }

    private void writeDrugData(final int patientId, final int treatmentId,
            @NotNull final BiopsyTreatmentDrugData drugData) {
        context.insertInto(DRUGDATA, DRUGDATA.STARTDATE, DRUGDATA.ENDDATE, DRUGDATA.NAME, DRUGDATA.TYPE,
                DRUGDATA.TREATMENTID, DRUGDATA.PATIENTID).values(Utils.toSQLDate(drugData.startDate()),
                Utils.toSQLDate(drugData.endDate()), drugData.name(), drugData.type(), treatmentId,
                patientId).execute();
    }

    private void writeTreatmentResponseData(final int patientId,
            @NotNull final BiopsyTreatmentResponseData treatmentResponseData) {
        context.insertInto(TREATMENTRESPONSEDATA, TREATMENTRESPONSEDATA.DATE, TREATMENTRESPONSEDATA.RESPONSE,
                TREATMENTRESPONSEDATA.MEASUREMENTDONE, TREATMENTRESPONSEDATA.TREATMENTID,
                TREATMENTRESPONSEDATA.PATIENTID).values(Utils.toSQLDate(treatmentResponseData.date()),
                treatmentResponseData.response(), treatmentResponseData.measurementDone(),
                treatmentResponseData.treatmentId(), patientId).execute();

    }

    private void writeSomaticVariantData(final int patientId, @NotNull final SomaticVariantData somaticVariantData) {
        context.insertInto(SOMATICVARIANTDATA, SOMATICVARIANTDATA.GENE, SOMATICVARIANTDATA.POSITION,
                SOMATICVARIANTDATA.REF, SOMATICVARIANTDATA.ALT, SOMATICVARIANTDATA.COSMICID,
                SOMATICVARIANTDATA.TOTALREADCOUNT, SOMATICVARIANTDATA.ALLELEREADCOUNT,
                SOMATICVARIANTDATA.PATIENTID).values(somaticVariantData.gene(), somaticVariantData.position(),
                somaticVariantData.ref(), somaticVariantData.alt(), somaticVariantData.cosmicID(),
                somaticVariantData.totalReadCount(), somaticVariantData.alleleReadCount(), patientId).execute();
    }
}
