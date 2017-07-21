package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.FORMSMETADATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSE;

import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentDrugData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.SampleData;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;
import org.jooq.Record;

class ClinicalDAO {

    @NotNull
    private final DSLContext context;

    ClinicalDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(PATIENT).execute();
        context.truncate(SAMPLE).execute();
        context.truncate(BIOPSY).execute();
        context.truncate(TREATMENT).execute();
        context.truncate(DRUG).execute();
        context.truncate(TREATMENTRESPONSE).execute();
        context.truncate(FORMSMETADATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    public void writeClinicalData(@NotNull final Patient patient) {
        final int patientId = writePatientData(patient.patientData());
        patient.sequencedBiopsies().forEach(biopsy -> writeSampleData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
    }

    private int writePatientData(@NotNull final PatientData patient) {
        final Record patientRecord = context.select(PATIENT.ID).from(PATIENT).where(PATIENT.CPCTID.eq(patient.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(PATIENT.ID);
        } else {
            final int patientId = context.insertInto(PATIENT, PATIENT.CPCTID, PATIENT.REGISTRATIONDATE, PATIENT.GENDER, PATIENT.ETHNICITY,
                    PATIENT.HOSPITAL, PATIENT.BIRTHYEAR, PATIENT.PRIMARYTUMORLOCATION, PATIENT.DEATHDATE)
                    .values(patient.cpctId(), Utils.toSQLDate(patient.registrationDate()), patient.gender(), patient.ethnicity(),
                            patient.hospital(), patient.birthYear(), patient.primaryTumorLocation(), Utils.toSQLDate(patient.deathDate()))
                    .returning(PATIENT.ID)
                    .fetchOne()
                    .getValue(PATIENT.ID);
            writeFormStatus(patientId, PATIENT.getName(), "demography", patient.demographyStatus(), patient.demographyLocked());
            writeFormStatus(patientId, PATIENT.getName(), "primaryTumor", patient.primaryTumorStatus(), patient.primaryTumorLocked());
            writeFormStatus(patientId, PATIENT.getName(), "eligibility", patient.eligibilityStatus(), patient.eligibilityLocked());
            writeFormStatus(patientId, PATIENT.getName(), "selectionCriteria", patient.selectionCriteriaStatus(),
                    patient.selectionCriteriaLocked());
            writeFormStatus(patientId, PATIENT.getName(), "death", patient.deathStatus(), patient.deathLocked());
            return patientId;
        }
    }

    private void writeSampleData(final int patientId, @NotNull final SampleData sample) {
        context.insertInto(SAMPLE, SAMPLE.SAMPLEID, SAMPLE.PATIENTID, SAMPLE.ARRIVALDATE)
                .values(sample.sampleId(), patientId, Utils.toSQLDate(sample.arrivalDate()))
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

    private void writeFormStatus(final int id, @NotNull final String tableName, @NotNull final String formName,
            @NotNull final String formStatus, @NotNull final String formLocked) {
        context.insertInto(FORMSMETADATA, FORMSMETADATA.ID, FORMSMETADATA.TABLENAME, FORMSMETADATA.FORM, FORMSMETADATA.STATUS,
                FORMSMETADATA.LOCKED).values(id, tableName, formName, formStatus, formLocked).execute();
    }
}
