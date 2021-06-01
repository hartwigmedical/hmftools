package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedDrug;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutablePreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.PreTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class PreTreatmentReader {

    private static final Logger LOGGER = LogManager.getLogger(PreTreatmentReader.class);

    private static final String STUDY_BASELINE = "SE.BASELINE";

    private static final String FORM_TREATMENT = "FRM.PRETHERAPY";
    private static final String ITEMGROUP_TREATMENT = "GRP.PRETHERAPY.PRETHERAPY";
    private static final String FIELD_PRETREATMENT_GIVEN = "FLD.PRETHERAPY.SYSTEMIC";
    private static final String FIELD_PRERADIOTHERAPY_GIVEN = "FLD.PRETHERAPY.RADIOTHER";

    private static final String ITEMGROUP_DRUGS = "GRP.PRETHERAPY.SYSTEMICTRT";
    private static final String FIELD_PRE_DRUG_START = "FLD.PRETHERAPY.SYSTEMICSTDTC";
    private static final String FIELD_PRE_DRUG_END = "FLD.PRETHERAPY.SYSTEMICENDTC";
    private static final String FIELD_PRE_DRUG = "FLD.PRETHERAPY.SYSTEMICREG";
    private static final String FIELD_PRE_BEST_RESPONSE = "FLD.PRETHERAPY.SYSTEMICRESP";

    @NotNull
    private final TreatmentCurator treatmentCurator;

    PreTreatmentReader(@NotNull final TreatmentCurator treatmentCurator) {
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    PreTreatmentData read(@NotNull EcrfPatient patient) {
        PreTreatmentData preTreatmentData = null;
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BASELINE)) {
            for (EcrfForm treatmentForm : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT)) {
                String treatmentGiven = readTreatmentGiven(treatmentForm);
                String radiotherapyGiven = readRadiotherapyGiven(treatmentForm);
                List<DrugData> drugs = readDrugs(treatmentForm);
                if (preTreatmentData == null) {
                    preTreatmentData = ImmutablePreTreatmentData.of(treatmentGiven, radiotherapyGiven, drugs, treatmentForm.status());
                } else {
                    LOGGER.warn("Multiple pre-therapy forms for found patient: {}", patient.patientId());
                }
            }
        }

        return preTreatmentData != null
                ? preTreatmentData
                : ImmutablePreTreatmentData.of(null, null, Lists.newArrayList(), FormStatus.undefined());
    }

    @NotNull
    private List<DrugData> readDrugs(@NotNull EcrfForm treatmentForm) {
        List<DrugData> drugs = Lists.newArrayList();
        for (EcrfItemGroup itemGroup : treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_DRUGS)) {
            LocalDate drugStart = itemGroup.readItemDate(FIELD_PRE_DRUG_START);
            LocalDate drugEnd = itemGroup.readItemDate(FIELD_PRE_DRUG_END);
            String drugName = itemGroup.readItemString(FIELD_PRE_DRUG);
            String bestResponse = itemGroup.readItemString(FIELD_PRE_BEST_RESPONSE);

            List<CuratedDrug> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
            drugs.add(ImmutableDrugData.of(drugName, drugStart, drugEnd, bestResponse, curatedDrugs));
        }
        return drugs;
    }

    @Nullable
    private static String readTreatmentGiven(@NotNull EcrfForm treatmentForm) {
        List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_PRETREATMENT_GIVEN);
        }
        return null;
    }

    @Nullable
    private static String readRadiotherapyGiven(@NotNull EcrfForm treatmentForm) {
        List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_PRERADIOTHERAPY_GIVEN);
        }
        return null;
    }
}
