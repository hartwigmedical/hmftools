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

import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PatientData;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;

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

    void writeClinicalData(@NotNull final Patient patient) {
        final int patientId = writePatientData(patient.patientData(), patient.preTreatmentData());
        patient.sequencedBiopsies().forEach(biopsy -> writeSampleData(patientId, biopsy));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
        patient.tumorMarkers().forEach(tumorMarker -> writeTumorMarkerData(patientId, tumorMarker));
    }

    private int writePatientData(@NotNull final PatientData patient, @NotNull PreTreatmentData preTreatmentData) {
        final Record patientRecord =
                context.select(PATIENT.ID).from(PATIENT).where(PATIENT.PATIENTIDENTIFIER.eq(patient.cpctId())).fetchOne();
        if (patientRecord != null) {
            return patientRecord.getValue(PATIENT.ID);
        } else {
            final int patientId = context.insertInto(PATIENT, PATIENT.PATIENTIDENTIFIER)
                    .values(patient.cpctId())
                    .returning(PATIENT.ID)
                    .fetchOne()
                    .getValue(PATIENT.ID);

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
                            patient.hospital(),
                            patient.birthYear(),
                            patient.cancerType().category(),
                            patient.cancerType().subcategory(),
                            Utils.toSQLDate(patient.deathDate()),
                            preTreatmentData.treatmentGiven(),
                            preTreatmentData.radiotherapyGiven(),
                            preTreatmentData.treatmentName())
                    .execute();
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "demography",
                    patient.demographyStatus().stateString(),
                    Boolean.toString(patient.demographyLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "primaryTumor",
                    patient.primaryTumorStatus().stateString(),
                    Boolean.toString(patient.primaryTumorLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "informedConsent",
                    patient.informedConsentStatus().stateString(),
                    Boolean.toString(patient.informedConsentLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "eligibility",
                    patient.eligibilityStatus().stateString(),
                    Boolean.toString(patient.eligibilityLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "selectionCriteria",
                    patient.selectionCriteriaStatus().stateString(),
                    Boolean.toString(patient.selectionCriteriaLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "death",
                    patient.deathStatus().stateString(),
                    Boolean.toString(patient.deathLocked()));
            writeFormStatus(patientId,
                    BASELINE.getName(),
                    "pretreatment",
                    preTreatmentData.formStatus().stateString(),
                    Boolean.toString(preTreatmentData.formLocked()));
            preTreatmentData.drugs()
                    .forEach(drug -> writePreTreatmentDrugData(patientId,
                            drug,
                            preTreatmentData.formStatus().stateString(),
                            Boolean.toString(preTreatmentData.formLocked())));
            return patientId;
        }
    }

    private void writePreTreatmentDrugData(final int patientId, @NotNull final DrugData drug, @NotNull final String formStatus,
            @NotNull final String formLocked) {
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
            writeFormStatus(id, PRETREATMENTDRUG.getName(), "pretreatment", formStatus, formLocked);
        });
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

    private void writeBiopsyData(final int patientId, @NotNull final BiopsyData biopsy) {
        context.insertInto(BIOPSY,
                BIOPSY.ID,
                BIOPSY.SAMPLEID,
                BIOPSY.PATIENTID,
                BIOPSY.BIOPSYTAKEN,
                BIOPSY.BIOPSYEVALUABLE,
                BIOPSY.BIOPSYSITE,
                BIOPSY.BIOPSYLOCATION,
                BIOPSY.BIOPSYDATE)
                .values(biopsy.id(),
                        biopsy.sampleId(),
                        patientId,
                        biopsy.biopsyTaken(),
                        biopsy.biopsyEvaluable(),
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
