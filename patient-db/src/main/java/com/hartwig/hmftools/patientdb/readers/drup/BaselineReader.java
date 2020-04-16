package com.hartwig.hmftools.patientdb.readers.drup;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.patientdb.curators.TumorLocationCurator;
import com.hartwig.hmftools.patientdb.data.BaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableBaselineData;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTumorLocation;

import org.jetbrains.annotations.NotNull;

class BaselineReader {

    private static final String STUDY_BASELINE = "SE.BAS";
    private static final String STUDY_REGISTRATION = "SE.REG";
    private static final String STUDY_END_OF_TRIAL = "SE.EOT";

    private static final String FORM_BASELINE = "FRM.BAS";
    private static final String FORM_CSF = "FRM.CSF"; // Not sure where CSF stands for, could be "clinical study fields"?
    private static final String FORM_REGISTRATION = "FRM.REG";
    private static final String FORM_END_OF_TRIAL = "FRM.EOT";

    private static final String ITEMGROUP_BASELINE = "GRP.BAS";
    private static final String ITEMGROUP_CSF = "GRP.CSF"; // Not sure where CSF stands for, could be "clinical study fields"?
    private static final String ITEMGROUP_REGISTRATION = "GRP.REG";
    private static final String ITEMGROUP_END_OF_TRIAL = "GRP.EOT";

    private static final String FIELD_INFORMED_CONSENT_DATE = "FLD.ICDTC";
    private static final String FIELD_HOSPITAL = "FLD.INST";
    private static final String FIELD_PRIMARY_TUMOR_LOCATION = "FLD.BASTTYP";
    private static final String FIELD_PRIMARY_TUMOR_LOCATION_OTHER = "FLD.BASTTOSP";
    private static final String FIELD_GENDER = "FLD.GEN";
    private static final String FIELD_BIRTH_YEAR = "FLD.YOB";
    private static final String FIELD_DEATH_DATE = "FLD.DEATHDTC";

    @NotNull
    private final TumorLocationCurator tumorLocationCurator;

    BaselineReader(@NotNull final TumorLocationCurator tumorLocationCurator) {
        this.tumorLocationCurator = tumorLocationCurator;
    }

    @NotNull
    BaselineData read(@NotNull EcrfPatient patient) {
        ImmutableBaselineData.Builder baselineBuilder = ImmutableBaselineData.builder()
                .curatedTumorLocation(ImmutableCuratedTumorLocation.of(null, null, null))
                .demographyStatus(FormStatus.undefined())
                .deathStatus(FormStatus.undefined())
                .eligibilityStatus(FormStatus.undefined())
                .selectionCriteriaStatus(FormStatus.undefined())
                .primaryTumorStatus(FormStatus.undefined())
                .informedConsentStatus(FormStatus.undefined());

        for (EcrfStudyEvent baselineEvent : patient.studyEventsPerOID(STUDY_BASELINE)) {
            setInformedConsentAndTumorLocation(baselineBuilder, baselineEvent);
        }

        for (EcrfStudyEvent registrationEvent : patient.studyEventsPerOID(STUDY_REGISTRATION)) {
            setBirthYearGenderHospital(baselineBuilder, registrationEvent);
        }

        for (EcrfStudyEvent endOfTrialEvent : patient.studyEventsPerOID(STUDY_END_OF_TRIAL)) {
            setDeathDate(baselineBuilder, endOfTrialEvent);
        }

        return baselineBuilder.build();
    }

    private void setInformedConsentAndTumorLocation(@NotNull ImmutableBaselineData.Builder builder, @NotNull EcrfStudyEvent studyEvent) {
        for (EcrfForm baselineForm : studyEvent.nonEmptyFormsPerOID(FORM_BASELINE)) {
            for (EcrfItemGroup baselineItemGroup : baselineForm.nonEmptyItemGroupsPerOID(ITEMGROUP_BASELINE)) {
                builder.informedConsentDate(baselineItemGroup.readItemDate(FIELD_INFORMED_CONSENT_DATE));

                String primaryTumorLocation = baselineItemGroup.readItemString(FIELD_PRIMARY_TUMOR_LOCATION);
                if (primaryTumorLocation != null && primaryTumorLocation.trim().toLowerCase().startsWith("other")) {
                    primaryTumorLocation = baselineItemGroup.readItemString(FIELD_PRIMARY_TUMOR_LOCATION_OTHER);
                }
                builder.curatedTumorLocation(tumorLocationCurator.search(primaryTumorLocation));
                // This is somewhat ugly, the states are too tied with CPCT datamodel.
                builder.informedConsentStatus(baselineForm.status());
                builder.primaryTumorStatus(baselineForm.status());
            }
        }
    }

    private void setBirthYearGenderHospital(@NotNull ImmutableBaselineData.Builder builder, @NotNull EcrfStudyEvent studyEvent) {
        for (EcrfForm csfForm : studyEvent.nonEmptyFormsPerOID(FORM_CSF)) {
            for (EcrfItemGroup csfItemGroup : csfForm.nonEmptyItemGroupsPerOID(ITEMGROUP_CSF)) {
                builder.gender(csfItemGroup.readItemString(FIELD_GENDER));

                String birthYear = csfItemGroup.readItemString(FIELD_BIRTH_YEAR);
                if (birthYear != null) {
                    builder.birthYear(Integer.parseInt(birthYear));
                }
                // This is somewhat ugly, the states are too tied with CPCT datamodel.
                builder.demographyStatus(csfForm.status());
                builder.eligibilityStatus(csfForm.status());
            }
        }

        for (EcrfForm registrationForm : studyEvent.nonEmptyFormsPerOID(FORM_REGISTRATION)) {
            for (EcrfItemGroup registrationItemGroup : registrationForm.nonEmptyItemGroupsPerOID(ITEMGROUP_REGISTRATION)) {
                builder.hospital(registrationItemGroup.readItemString(FIELD_HOSPITAL));
                // This is somewhat ugly, the states are too tied with CPCT datamodel.
                builder.selectionCriteriaStatus(registrationForm.status());
            }
        }
    }

    private void setDeathDate(@NotNull ImmutableBaselineData.Builder builder, @NotNull EcrfStudyEvent studyEvent) {
        for (EcrfForm endOfTrialForm : studyEvent.nonEmptyFormsPerOID(FORM_END_OF_TRIAL)) {
            for (EcrfItemGroup endOfTrialItemGroup : endOfTrialForm.nonEmptyItemGroupsPerOID(ITEMGROUP_END_OF_TRIAL)) {
                builder.deathDate(endOfTrialItemGroup.readItemDate(FIELD_DEATH_DATE));
                builder.deathStatus(endOfTrialForm.status());
            }
        }
    }
}
