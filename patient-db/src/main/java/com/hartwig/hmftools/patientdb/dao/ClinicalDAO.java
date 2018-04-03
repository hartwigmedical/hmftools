package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.FORMSMETADATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PRETREATMENTDRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TUMORMARKER;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;

import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class ClinicalDAO {

    @NotNull
    private final DSLContext context;

    ClinicalDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(PATIENT).execute();
        context.truncate(BASELINE).execute();
        context.truncate(PRETREATMENTDRUG).execute();
        context.truncate(SAMPLE).execute();
        context.truncate(BIOPSY).execute();
        context.truncate(TREATMENT).execute();
        context.truncate(DRUG).execute();
        context.truncate(TREATMENTRESPONSE).execute();
        context.truncate(TUMORMARKER).execute();
        context.truncate(FORMSMETADATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writeFullClinicalData(@NotNull Patient patient) {
        final int patientId = writePatientIdentifier(patient.patientIdentifier());
        writeBaselineData(patientId, patient.baselineData(), patient.preTreatmentData());
        patient.sequencedBiopsies().forEach(sample -> writeSampleData(patientId, sample));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
        patient.tumorMarkers().forEach(tumorMarker -> writeTumorMarkerData(patientId, tumorMarker));
    }

    void writeSampleClinicalData(@NotNull String patientIdentifier, @NotNull List<SampleData> samples) {
        final int patientId = writePatientIdentifier(patientIdentifier);
        samples.forEach(sample -> writeSampleData(patientId, sample));
    }

    private int writePatientIdentifier(@NotNull String patientIdentifier) {
        return context.insertInto(PATIENT, PATIENT.PATIENTIDENTIFIER)
                .values(patientIdentifier)
                .returning(PATIENT.ID)
                .fetchOne()
                .getValue(PATIENT.ID);
    }

    private void writeSampleData(final int patientId, @NotNull final SampleData sample) {
        context.insertInto(SAMPLE,
                SAMPLE.SAMPLEID,
                SAMPLE.PATIENTID,
                SAMPLE.ARRIVALDATE,
                SAMPLE.SAMPLINGDATE,
                SAMPLE.DNANANOGRAMS,
                SAMPLE.TUMORPERCENTAGE)
                .values(sample.sampleId(),
                        patientId,
                        Utils.toSQLDate(sample.arrivalDate()),
                        Utils.toSQLDate(sample.samplingDate()),
                        sample.dnaNanograms(),
                        sample.tumorPercentage())
                .execute();
    }

    private void writeBaselineData(int patientId, @NotNull BaselineData patient,
            @NotNull PreTreatmentData preTreatmentData) {
        context.insertInto(BASELINE,
                BASELINE.PATIENTID,
                BASELINE.REGISTRATIONDATE,
                BASELINE.INFORMEDCONSENTDATE,
                BASELINE.GENDER,
                BASELINE.HOSPITAL,
                BASELINE.BIRTHYEAR,
                BASELINE.CANCERTYPE,
                BASELINE.CANCERSUBTYPE,
                BASELINE.DEATHDATE,
                BASELINE.HASSYSTEMICPRETREATMENT,
                BASELINE.HASRADIOTHERAPYPRETREATMENT,
                BASELINE.PRETREATMENTS)
                .values(patientId,
                        Utils.toSQLDate(patient.registrationDate()),
                        Utils.toSQLDate(patient.informedConsentDate()),
                        patient.gender(),
                        patient.hospital(), patient.birthYear(), patient.cancerType().type(), patient.cancerType().subType(),
                        Utils.toSQLDate(patient.deathDate()),
                        preTreatmentData.treatmentGiven(),
                        preTreatmentData.radiotherapyGiven(),
                        preTreatmentData.treatmentName())
                .execute();

        preTreatmentData.drugs()
                .forEach(drug -> writePreTreatmentDrugData(patientId, drug, preTreatmentData.formStatus(), preTreatmentData.formLocked()));

        writeBaselineFormStatus(patientId, "demography", patient.demographyStatus(), patient.demographyLocked());
        writeBaselineFormStatus(patientId, "primaryTumor", patient.primaryTumorStatus(), patient.primaryTumorLocked());
        writeBaselineFormStatus(patientId, "informedConsent", patient.informedConsentStatus(), patient.informedConsentLocked());
        writeBaselineFormStatus(patientId, "eligibility", patient.eligibilityStatus(), patient.eligibilityLocked());
        writeBaselineFormStatus(patientId, "selectionCriteria", patient.selectionCriteriaStatus(), patient.selectionCriteriaLocked());
        writeBaselineFormStatus(patientId, "death", patient.deathStatus(), patient.deathLocked());
        writeBaselineFormStatus(patientId, "pretreatment", preTreatmentData.formStatus(), preTreatmentData.formLocked());
    }

    private void writeBaselineFormStatus(int patientId, @NotNull String form, @NotNull FormStatusState formStatusState, boolean locked) {
        writeFormStatus(patientId, BASELINE.getName(), form, formStatusState.stateString(), Boolean.toString(locked));
    }

    private void writePreTreatmentDrugData(int patientId, @NotNull DrugData drug, @NotNull FormStatusState formStatus, boolean formLocked) {
        drug.filteredCuratedTreatments().forEach(curatedTreatment -> {
            final int id = context.insertInto(PRETREATMENTDRUG,
                    PRETREATMENTDRUG.PATIENTID,
                    PRETREATMENTDRUG.STARTDATE,
                    PRETREATMENTDRUG.ENDDATE,
                    PRETREATMENTDRUG.NAME,
                    PRETREATMENTDRUG.TYPE,
                    PRETREATMENTDRUG.BESTRESPONSE)
                    .values(patientId,
                            Utils.toSQLDate(drug.startDate()),
                            Utils.toSQLDate(drug.endDate()),
                            curatedTreatment.name(),
                            curatedTreatment.type(),
                            drug.bestResponse())
                    .returning(PRETREATMENTDRUG.ID)
                    .fetchOne()
                    .getValue(PRETREATMENTDRUG.ID);

            writeFormStatus(id, PRETREATMENTDRUG.getName(), "pretreatment", formStatus.stateString(), Boolean.toString(formLocked));
        });
    }

    private void writeBiopsyData(final int patientId, @NotNull final BiopsyData biopsy) {
        context.insertInto(BIOPSY,
                BIOPSY.ID,
                BIOPSY.SAMPLEID,
                BIOPSY.PATIENTID,
                BIOPSY.BIOPSYTAKEN,
                BIOPSY.BIOPSYEVALUABLE,
                BIOPSY.BIOPSYTYPE,
                BIOPSY.BIOPSYSITE,
                BIOPSY.BIOPSYLOCATION,
                BIOPSY.BIOPSYDATE)
                .values(biopsy.id(),
                        biopsy.sampleId(),
                        patientId,
                        biopsy.biopsyTaken(),
                        biopsy.biopsyEvaluable(),
                        biopsy.curatedType(),
                        biopsy.site(),
                        biopsy.location(),
                        Utils.toSQLDate(biopsy.date()))
                .execute();
        writeFormStatus(biopsy.id(), BIOPSY.getName(), "biopsy", biopsy.formStatus().stateString(), Boolean.toString(biopsy.formLocked()));
    }

    private void writeTreatmentData(final int patientId, @NotNull final BiopsyTreatmentData treatment) {
        context.insertInto(TREATMENT,
                TREATMENT.ID,
                TREATMENT.BIOPSYID,
                TREATMENT.PATIENTID,
                TREATMENT.TREATMENTGIVEN,
                TREATMENT.RADIOTHERAPYGIVEN,
                TREATMENT.STARTDATE,
                TREATMENT.ENDDATE,
                TREATMENT.NAME,
                TREATMENT.TYPE)
                .values(treatment.id(),
                        treatment.biopsyId(),
                        patientId,
                        treatment.treatmentGiven(),
                        treatment.radiotherapyGiven(),
                        Utils.toSQLDate(treatment.startDate()),
                        Utils.toSQLDate(treatment.endDate()),
                        treatment.treatmentName(),
                        treatment.type())
                .execute();
        writeFormStatus(treatment.id(),
                TREATMENT.getName(),
                "treatment",
                treatment.formStatus().stateString(),
                Boolean.toString(treatment.formLocked()));
        treatment.drugs()
                .forEach(drug -> writeDrugData(patientId,
                        treatment.id(),
                        drug,
                        treatment.formStatus().stateString(),
                        Boolean.toString(treatment.formLocked())));
    }

    private void writeDrugData(final int patientId, final int treatmentId, @NotNull final DrugData drug, @NotNull final String formStatus,
            @NotNull final String formLocked) {
        drug.filteredCuratedTreatments().forEach(curatedTreatment -> {
            final int id = context.insertInto(DRUG, DRUG.TREATMENTID, DRUG.PATIENTID, DRUG.STARTDATE, DRUG.ENDDATE, DRUG.NAME, DRUG.TYPE)
                    .values(treatmentId,
                            patientId,
                            Utils.toSQLDate(drug.startDate()),
                            Utils.toSQLDate(drug.endDate()),
                            curatedTreatment.name(),
                            curatedTreatment.type())
                    .returning(DRUG.ID)
                    .fetchOne()
                    .getValue(DRUG.ID);
            writeFormStatus(id, DRUG.getName(), "treatment", formStatus, formLocked);
        });
    }

    private void writeTreatmentResponseData(final int patientId, @NotNull final BiopsyTreatmentResponseData treatmentResponse) {
        final int id = context.insertInto(TREATMENTRESPONSE,
                TREATMENTRESPONSE.TREATMENTID,
                TREATMENTRESPONSE.PATIENTID,
                TREATMENTRESPONSE.RESPONSEDATE,
                TREATMENTRESPONSE.RESPONSE,
                TREATMENTRESPONSE.MEASUREMENTDONE,
                TREATMENTRESPONSE.BONEONLYDISEASE)
                .values(treatmentResponse.treatmentId(),
                        patientId,
                        Utils.toSQLDate(treatmentResponse.date()),
                        treatmentResponse.response(),
                        treatmentResponse.measurementDone(),
                        treatmentResponse.boneOnlyDisease())
                .returning(TREATMENTRESPONSE.ID)
                .fetchOne()
                .getValue(TREATMENTRESPONSE.ID);

        writeFormStatus(id,
                TREATMENTRESPONSE.getName(),
                "treatmentResponse",
                treatmentResponse.formStatus().stateString(),
                Boolean.toString(treatmentResponse.formLocked()));
    }

    private void writeTumorMarkerData(final int patientId, @NotNull final TumorMarkerData tumorMarker) {
        final int id = context.insertInto(TUMORMARKER,
                TUMORMARKER.PATIENTID,
                TUMORMARKER.DATE,
                TUMORMARKER.MARKER,
                TUMORMARKER.MEASUREMENT,
                TUMORMARKER.UNIT)
                .values(patientId, Utils.toSQLDate(tumorMarker.date()), tumorMarker.marker(), tumorMarker.measurement(), tumorMarker.unit())
                .returning(TUMORMARKER.ID)
                .fetchOne()
                .getValue(TUMORMARKER.ID);

        writeFormStatus(id,
                TUMORMARKER.getName(),
                "tumorMarker",
                tumorMarker.formStatus().stateString(),
                Boolean.toString(tumorMarker.formLocked()));
    }

    private void writeFormStatus(final int id, @NotNull final String tableName, @NotNull final String formName,
            @NotNull final String formStatus, @NotNull final String formLocked) {
        context.insertInto(FORMSMETADATA,
                FORMSMETADATA.ID,
                FORMSMETADATA.TABLENAME,
                FORMSMETADATA.FORM,
                FORMSMETADATA.STATUS,
                FORMSMETADATA.LOCKED).values(id, tableName, formName, formStatus, formLocked).execute();
    }
}
