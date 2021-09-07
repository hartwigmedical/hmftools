package com.hartwig.hmftools.patientdb.clinical.readers.cpct;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.patientdb.clinical.curators.TreatmentCurator;
import com.hartwig.hmftools.patientdb.clinical.datamodel.BiopsyTreatmentData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedDrug;
import com.hartwig.hmftools.patientdb.clinical.datamodel.DrugData;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableDrugData;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfStudyEvent;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class BiopsyTreatmentReader {

    private static final String STUDY_AFTERBIOPT = "SE.AFTERBIOPT";
    private static final String FORM_TREATMENT = "FRM.TRTAFTER";

    private static final String ITEMGROUP_TREATMENT_AFTER = "GRP.TRTAFTER.TRTAFTER";
    private static final String FIELD_TREATMENT_GIVEN = "FLD.TRTAFTER.SYSTEMICST";
    private static final String FIELD_RADIOTHERAPY_GIVEN = "FLD.TRTAFTER.RADIOTHERST";

    private static final String ITEMGROUP_SYSPOSTBIO = "GRP.TRTAFTER.SYSPOSTBIO";
    private static final String FIELD_DRUG_START = "FLD.TRTAFTER.SYSSTDT";
    private static final String FIELD_DRUG_END = "FLD.TRTAFTER.SYSENDT";
    private static final String FIELD_DRUG = "FLD.TRTAFTER.PLANNEDTRT";
    private static final String FIELD_DRUG_OTHER = "FLD.TRTAFTER.SYSREGPOST";

    @NotNull
    private final TreatmentCurator treatmentCurator;

    BiopsyTreatmentReader(@NotNull final TreatmentCurator treatmentCurator) {
        this.treatmentCurator = treatmentCurator;
    }

    @NotNull
    List<BiopsyTreatmentData> read(@NotNull EcrfPatient patient) {
        List<BiopsyTreatmentData> treatments = Lists.newArrayList();
        for (EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_AFTERBIOPT)) {
            for (EcrfForm treatmentForm : studyEvent.nonEmptyFormsPerOID(FORM_TREATMENT)) {
                String treatmentGiven = readTreatmentGiven(treatmentForm);
                String radiotherapyGiven = readRadiotherapyGiven(treatmentForm);
                List<DrugData> drugs = readDrugs(treatmentForm);
                treatments.add(BiopsyTreatmentData.of(null, treatmentGiven, radiotherapyGiven, drugs, treatmentForm.status()));
            }
        }
        return treatments;
    }

    @NotNull
    private List<DrugData> readDrugs(@NotNull EcrfForm treatmentForm) {
        List<DrugData> drugs = Lists.newArrayList();
        for (EcrfItemGroup itemGroup : treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_SYSPOSTBIO)) {
            LocalDate drugStart = itemGroup.readItemDate(FIELD_DRUG_START);
            LocalDate drugEnd = itemGroup.readItemDate(FIELD_DRUG_END);
            String drugName = itemGroup.readItemString(FIELD_DRUG);
            if (drugName == null || drugName.trim().toLowerCase().startsWith("other")) {
                drugName = itemGroup.readItemString(FIELD_DRUG_OTHER);
            }
            if (drugName != null || drugStart != null || drugEnd != null) {
                List<CuratedDrug> curatedDrugs = drugName == null ? Lists.newArrayList() : treatmentCurator.search(drugName);
                drugs.add(ImmutableDrugData.of(drugName, drugStart, drugEnd, null, curatedDrugs));
            }
        }
        return drugs;
    }

    @Nullable
    private static String readTreatmentGiven(@NotNull EcrfForm treatmentForm) {
        List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT_AFTER);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_TREATMENT_GIVEN);
        }
        return null;
    }

    @Nullable
    private static String readRadiotherapyGiven(@NotNull EcrfForm treatmentForm) {
        List<EcrfItemGroup> itemGroups = treatmentForm.nonEmptyItemGroupsPerOID(ITEMGROUP_TREATMENT_AFTER);
        if (itemGroups.size() > 0) {
            return itemGroups.get(0).readItemString(FIELD_RADIOTHERAPY_GIVEN);
        }
        return null;
    }
}
