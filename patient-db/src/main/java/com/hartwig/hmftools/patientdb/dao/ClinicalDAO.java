package com.hartwig.hmftools.patientdb.dao;

import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BASELINE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.BIOPSY;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DOIDNODE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.DRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.FORMSMETADATA;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PATIENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.PRETREATMENTDRUG;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.RANOMEASUREMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SAMPLE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.SNOMED;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENT;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TREATMENTRESPONSE;
import static com.hartwig.hmftools.patientdb.database.hmfpatients.Tables.TUMORMARKER;

import java.util.List;

import com.hartwig.hmftools.common.doid.DoidNode;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.BiopsyData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.data.BiopsyTreatmentResponseData;
import com.hartwig.hmftools.patientdb.data.DrugData;
import com.hartwig.hmftools.patientdb.data.Patient;
import com.hartwig.hmftools.patientdb.data.PreTreatmentData;
import com.hartwig.hmftools.patientdb.data.RanoMeasurementData;
import com.hartwig.hmftools.patientdb.data.SampleData;
import com.hartwig.hmftools.patientdb.data.TumorMarkerData;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jooq.DSLContext;

class ClinicalDAO {

    private static final Logger LOGGER = LogManager.getLogger(ClinicalDAO.class);

    @NotNull
    private final DSLContext context;

    ClinicalDAO(@NotNull final DSLContext context) {
        this.context = context;
    }

    void clear() {
        context.execute("SET FOREIGN_KEY_CHECKS = 0;");
        context.truncate(PATIENT).execute();
        context.truncate(BASELINE).execute();
        context.truncate(DOIDNODE).execute();
        context.truncate(SNOMED).execute();
        context.truncate(PRETREATMENTDRUG).execute();
        context.truncate(SAMPLE).execute();
        context.truncate(BIOPSY).execute();
        context.truncate(TREATMENT).execute();
        context.truncate(DRUG).execute();
        context.truncate(TREATMENTRESPONSE).execute();
        context.truncate(TUMORMARKER).execute();
        context.truncate(RANOMEASUREMENT).execute();
        context.truncate(FORMSMETADATA).execute();
        context.execute("SET FOREIGN_KEY_CHECKS = 1;");
    }

    void writeFullClinicalData(@NotNull Patient patient, boolean blacklisted) {
        int patientId = writePatient(patient.patientIdentifier(), blacklisted);
        writePatientData(patientId, patient.baselineData(), patient.preTreatmentData());
        patient.sequencedBiopsies().forEach(sample -> writeSampleData(patientId, sample));
        patient.clinicalBiopsies().forEach(biopsy -> writeBiopsyData(patientId, biopsy));
        patient.treatments().forEach(treatment -> writeTreatmentData(patientId, treatment));
        patient.treatmentResponses().forEach(response -> writeTreatmentResponseData(patientId, response));
        patient.tumorMarkers().forEach(tumorMarker -> writeTumorMarkerData(patientId, tumorMarker));
        patient.ranoMeasurements().forEach(ranoMeasurement -> writeRanoMeasurementData(patientId, ranoMeasurement));
    }

    void writeSampleClinicalData(@NotNull String patientIdentifier, boolean blacklisted, @NotNull List<SampleData> samples) {
        int patientId = writePatient(patientIdentifier, blacklisted);
        samples.forEach(sample -> writeSampleData(patientId, sample));
    }

    private int writePatient(@NotNull String patientIdentifier, boolean blacklisted) {
        return context.insertInto(PATIENT, PATIENT.PATIENTIDENTIFIER, PATIENT.BLACKLISTED)
                .values(patientIdentifier, (byte) (blacklisted ? 1 : 0))
                .returning(PATIENT.ID)
                .fetchOne()
                .getValue(PATIENT.ID);
    }

    private void writeSampleData(int patientId, @NotNull SampleData sample) {
        context.insertInto(SAMPLE,
                SAMPLE.SAMPLEID,
                SAMPLE.PATIENTID,
                SAMPLE.SETNAME,
                SAMPLE.ARRIVALDATE,
                SAMPLE.SAMPLINGDATE,
                SAMPLE.DNANANOGRAMS,
                SAMPLE.LIMSPRIMARYTUMOR,
                SAMPLE.PATHOLOGYTUMORPERCENTAGE)
                .values(sample.sampleId(),
                        patientId,
                        sample.setName(),
                        Utils.toSQLDate(sample.arrivalDate()),
                        Utils.toSQLDate(sample.samplingDate()),
                        sample.dnaNanograms(),
                        sample.limsPrimaryTumor(),
                        sample.pathologyTumorPercentage())
                .execute();
    }

    private void writePatientData(int patientId, @NotNull BaselineData patient, @NotNull PreTreatmentData preTreatmentData) {
        // preTreatmentTypes exceeds the usual 255 length of varchar fields in production.
        String preTreatmentTypes = preTreatmentData.concatenatedType();
        if (preTreatmentTypes != null && preTreatmentTypes.length() > BASELINE.PRETREATMENTSTYPE.getDataType().length()) {
            LOGGER.warn("Truncating pre-treatment type: {}", preTreatmentTypes);
            preTreatmentTypes = preTreatmentTypes.substring(0, BASELINE.PRETREATMENTSTYPE.getDataType().length());
        }

        String preTreatmentMechanism = preTreatmentData.concatenatedMechanism();
        if (preTreatmentMechanism != null && preTreatmentMechanism.length() > BASELINE.PRETREATMENTSMECHANISM.getDataType().length()) {
            LOGGER.warn("Truncating pre-treatment type: {}", preTreatmentMechanism);
            preTreatmentMechanism = preTreatmentMechanism.substring(0, BASELINE.PRETREATMENTSMECHANISM.getDataType().length());
        }
        byte primaryTumorOverridden = 0;

        context.insertInto(BASELINE,
                BASELINE.PATIENTID,
                BASELINE.REGISTRATIONDATE,
                BASELINE.INFORMEDCONSENTDATE,
                BASELINE.GENDER,
                BASELINE.HOSPITAL,
                BASELINE.BIRTHYEAR,
                BASELINE.PRIMARYTUMORLOCATION,
                BASELINE.PRIMARYTUMORSUBLOCATION,
                BASELINE.PRIMARYTUMORTYPE,
                BASELINE.PRIMARYTUMORSUBTYPE,
                BASELINE.PRIMARYTUMOREXTRADETAILS,
                BASELINE.PRIMARYTUMOROVERRIDDEN,
                BASELINE.DEATHDATE,
                BASELINE.HASSYSTEMICPRETREATMENT,
                BASELINE.HASRADIOTHERAPYPRETREATMENT,
                BASELINE.PRETREATMENTS,
                BASELINE.PRETREATMENTSTYPE,
                BASELINE.PRETREATMENTSMECHANISM)
                .values(patientId,
                        Utils.toSQLDate(patient.registrationDate()),
                        Utils.toSQLDate(patient.informedConsentDate()),
                        patient.gender(),
                        patient.hospital(),
                        patient.birthYear(),
                        patient.curatedPrimaryTumor().location(),
                        patient.curatedPrimaryTumor().subLocation(),
                        patient.curatedPrimaryTumor().type(),
                        patient.curatedPrimaryTumor().subType(),
                        patient.curatedPrimaryTumor().extraDetails(),
                        primaryTumorOverridden,
                        Utils.toSQLDate(patient.deathDate()),
                        preTreatmentData.treatmentGiven(),
                        preTreatmentData.radiotherapyGiven(),
                        preTreatmentData.treatmentName(),
                        preTreatmentTypes,
                        preTreatmentMechanism)
                .execute();

        List<DoidNode> doidNodes = patient.curatedPrimaryTumor().doidNodes();
        if (doidNodes != null) {
            for (DoidNode doidNode : doidNodes) {
                writeDoidNode(patientId, doidNode);
            }
        }

        List<String> snomedConceptIds = patient.curatedPrimaryTumor().snomedConceptIds();
        if (snomedConceptIds != null) {
            for (String snomedConceptId : snomedConceptIds) {
                writeSnomedConceptId(patientId, snomedConceptId);
            }
        }

        for (DrugData drug : preTreatmentData.drugs()) {
            writePreTreatmentDrugData(patientId, drug, preTreatmentData.formStatus());
        }

        writeBaselineFormStatus(patientId, "demography", patient.demographyStatus());
        writeBaselineFormStatus(patientId, "primaryTumor", patient.primaryTumorStatus());
        writeBaselineFormStatus(patientId, "informedConsent", patient.informedConsentStatus());
        writeBaselineFormStatus(patientId, "eligibility", patient.eligibilityStatus());
        writeBaselineFormStatus(patientId, "selectionCriteria", patient.selectionCriteriaStatus());
        writeBaselineFormStatus(patientId, "death", patient.deathStatus());
        writeBaselineFormStatus(patientId, "pretreatment", preTreatmentData.formStatus());
    }

    private void writeBaselineFormStatus(int patientId, @NotNull String form, @NotNull FormStatus formStatus) {
        writeFormStatus(patientId, BASELINE.getName(), form, formStatus);
    }

    private void writeDoidNode(int patientId, @NotNull DoidNode doidNode) {
        context.insertInto(DOIDNODE, DOIDNODE.PATIENTID, DOIDNODE.DOID, DOIDNODE.DOIDTERM, DOIDNODE.SNOMEDCONCEPTID)
                .values(patientId, doidNode.doid(), doidNode.doidTerm(), doidNode.snomedConceptId())
                .execute();
    }

    private void writeSnomedConceptId(int patientId, @NotNull String snomedConceptId) {
        context.insertInto(SNOMED, SNOMED.PATIENTID, SNOMED.SNOMEDCONCEPTID).values(patientId, snomedConceptId).execute();
    }

    private void writePreTreatmentDrugData(int patientId, @NotNull DrugData drug, @NotNull FormStatus formStatus) {
        drug.filteredCuratedDrugs().forEach(curatedTreatment -> {
            int id = context.insertInto(PRETREATMENTDRUG,
                    PRETREATMENTDRUG.PATIENTID,
                    PRETREATMENTDRUG.STARTDATE,
                    PRETREATMENTDRUG.ENDDATE,
                    PRETREATMENTDRUG.NAME,
                    PRETREATMENTDRUG.TYPE,
                    PRETREATMENTDRUG.MECHANISM,
                    PRETREATMENTDRUG.BESTRESPONSE)
                    .values(patientId,
                            Utils.toSQLDate(drug.startDate()),
                            Utils.toSQLDate(drug.endDate()),
                            curatedTreatment.name(),
                            curatedTreatment.type(),
                            curatedTreatment.mechanism(),
                            drug.bestResponse())
                    .returning(PRETREATMENTDRUG.ID)
                    .fetchOne()
                    .getValue(PRETREATMENTDRUG.ID);

            writeFormStatus(id, PRETREATMENTDRUG.getName(), "pretreatment", formStatus);
        });
    }

    private void writeBiopsyData(int patientId, @NotNull BiopsyData biopsy) {
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
        writeFormStatus(biopsy.id(), BIOPSY.getName(), "biopsy", biopsy.formStatus());
    }

    private void writeTreatmentData(int patientId, @NotNull BiopsyTreatmentData treatment) {
        context.insertInto(TREATMENT,
                TREATMENT.ID,
                TREATMENT.BIOPSYID,
                TREATMENT.PATIENTID,
                TREATMENT.TREATMENTGIVEN,
                TREATMENT.RADIOTHERAPYGIVEN,
                TREATMENT.STARTDATE,
                TREATMENT.ENDDATE,
                TREATMENT.NAME,
                TREATMENT.TYPE,
                TREATMENT.MECHANISM)
                .values(treatment.id(),
                        treatment.biopsyId(),
                        patientId,
                        treatment.treatmentGiven(),
                        treatment.radiotherapyGiven(),
                        Utils.toSQLDate(treatment.startDate()),
                        Utils.toSQLDate(treatment.endDate()),
                        treatment.treatmentName(),
                        treatment.consolidatedType(),
                        treatment.consolidatedMechanism())
                .execute();
        writeFormStatus(treatment.id(), TREATMENT.getName(), "treatment", treatment.formStatus());
        treatment.drugs().forEach(drug -> writeDrugData(patientId, treatment.id(), drug, treatment.formStatus()));
    }

    private void writeDrugData(int patientId, int treatmentId, @NotNull DrugData drug, @NotNull FormStatus formStatus) {
        drug.filteredCuratedDrugs().forEach(curatedTreatment -> {
            int id = context.insertInto(DRUG,
                    DRUG.TREATMENTID,
                    DRUG.PATIENTID,
                    DRUG.STARTDATE,
                    DRUG.ENDDATE,
                    DRUG.NAME,
                    DRUG.TYPE,
                    DRUG.MECHANISM)
                    .values(treatmentId,
                            patientId,
                            Utils.toSQLDate(drug.startDate()),
                            Utils.toSQLDate(drug.endDate()),
                            curatedTreatment.name(),
                            curatedTreatment.type(),
                            curatedTreatment.mechanism())
                    .returning(DRUG.ID)
                    .fetchOne()
                    .getValue(DRUG.ID);
            writeFormStatus(id, DRUG.getName(), "treatment", formStatus);
        });
    }

    private void writeTreatmentResponseData(int patientId, @NotNull BiopsyTreatmentResponseData treatmentResponse) {
        int id = context.insertInto(TREATMENTRESPONSE,
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

        writeFormStatus(id, TREATMENTRESPONSE.getName(), "treatmentResponse", treatmentResponse.formStatus());
    }

    private void writeTumorMarkerData(int patientId, @NotNull TumorMarkerData tumorMarker) {
        int id = context.insertInto(TUMORMARKER,
                TUMORMARKER.PATIENTID,
                TUMORMARKER.DATE,
                TUMORMARKER.MARKER,
                TUMORMARKER.MEASUREMENT,
                TUMORMARKER.UNIT)
                .values(patientId, Utils.toSQLDate(tumorMarker.date()), tumorMarker.marker(), tumorMarker.measurement(), tumorMarker.unit())
                .returning(TUMORMARKER.ID)
                .fetchOne()
                .getValue(TUMORMARKER.ID);

        writeFormStatus(id, TUMORMARKER.getName(), "tumorMarker", tumorMarker.formStatus());
    }

    private void writeRanoMeasurementData(int patientId, @NotNull RanoMeasurementData RanoMeasurement) {
        int id = context.insertInto(RANOMEASUREMENT,
                RANOMEASUREMENT.PATIENTID,
                RANOMEASUREMENT.RESPONSEDATE,
                RANOMEASUREMENT.THERAPYGIVEN,
                RANOMEASUREMENT.TARGETLESIONRESPONSE,
                RANOMEASUREMENT.NOTARGETLESIONRESPONSE,
                RANOMEASUREMENT.OVERALLRESPONSE)
                .values(patientId,
                        Utils.toSQLDate(RanoMeasurement.responseDate()),
                        RanoMeasurement.therapyGiven(),
                        RanoMeasurement.targetLesionResponse(),
                        RanoMeasurement.noTargetLesionResponse(),
                        RanoMeasurement.overallResponse())
                .returning(RANOMEASUREMENT.ID)
                .fetchOne()
                .getValue(RANOMEASUREMENT.ID);

        writeFormStatus(id, RANOMEASUREMENT.getName(), "RanoMeasurement", RanoMeasurement.formStatus());
    }

    private void writeFormStatus(int id, @NotNull String tableName, @NotNull String formName, @NotNull FormStatus formStatus) {
        context.insertInto(FORMSMETADATA,
                FORMSMETADATA.ID,
                FORMSMETADATA.TABLENAME,
                FORMSMETADATA.FORM,
                FORMSMETADATA.STATUS,
                FORMSMETADATA.LOCKED).values(id, tableName, formName, formStatus.stateString(), formStatus.lockedString()).execute();
    }
}
