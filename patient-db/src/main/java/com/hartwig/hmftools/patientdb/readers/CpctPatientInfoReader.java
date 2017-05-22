package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.Map;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfForm;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfItemGroup;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfStudyEvent;
import com.hartwig.hmftools.patientdb.data.PatientInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CpctPatientInfoReader {
    private static final Logger LOGGER = LogManager.getLogger(CpctPatientInfoReader.class);

    private static final String STUDY_BASELINE = "SE.BASELINE";
    private static final String STUDY_ENDSTUDY = "SE.ENDSTUDY";

    private static final String FORM_DEMOGRAPHY = "FRM.DEMOGRAPHY";
    private static final String FORM_CARCINOMA = "FRM.CARCINOMA";
    private static final String FORM_ELIGIBILITY = "FRM.ELIGIBILITY";
    private static final String FORM_SELCRIT = "FRM.SELCRIT";
    private static final String FORM_DEATH = "FRM.DEATH";

    private static final String ITEMGROUP_DEMOGRAPHY = "GRP.DEMOGRAPHY.DEMOGRAPHY";
    private static final String ITEMGROUP_CARCINOMA = "GRP.CARCINOMA.CARCINOMA";
    private static final String ITEMGROUP_ELIGIBILITY = "GRP.ELIGIBILITY.ELIGIBILITY";
    private static final String ITEMGROUP_SELCRIT = "GRP.SELCRIT.SELCRIT";
    private static final String ITEMGROUP_DEATH = "GRP.DEATH.DEATH";

    private static final String FIELD_SEX = "FLD.DEMOGRAPHY.SEX";
    private static final String FIELD_ETHNICITY = "FLD.DEMOGRAPHY.ETHNIC";
    //private static final String FIELD_ETHNICITY_SPECIFY = "FLD.DEMOGRAPHY.ETHNICSP";

    private static final String FIELD_REGISTRATION_DATE = "FLD.ELIGIBILITY.REGDTC";
    private static final String FIELD_HOSPITAL1 = "FLD.ELIGIBILITY.HOSPITAL";
    private static final String FIELD_BIRTH_YEAR2 = "FLD.ELIGIBILITY.BIRTHYEAR";
    private static final String FIELD_BIRTH_YEAR3 = "FLD.ELIGIBILITY.BIRTHDTCES";

    private static final String FIELD_HOSPITAL2 = "FLD.SELCRIT.NHOSPITAL";
    private static final String FIELD_BIRTH_YEAR1 = "FLD.SELCRIT.NBIRTHYEAR";

    private static final String FIELD_TUMOR_LOCATION = "FLD.CARCINOMA.PTUMLOC";

    private static final String FIELD_DEATH_DATE = "FLD.DEATH.DDEATHDTC";

    private static final String DATAMODEL_HOSPITAL1 = "BASELINE.ELIGIBILITY.ELIGIBILITY.HOSPITAL";

    private static final String DATAMODEL_HOSPITAL2 = "BASELINE.SELCRIT.SELCRIT.NHOSPITAL";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final Map<Integer, String> hospitals;

    public CpctPatientInfoReader(@NotNull final CpctEcrfModel model) {
        this.hospitals = getHospitals(model);
    }

    @NotNull
    public PatientInfo read(@NotNull final EcrfPatient patient) {
        LOGGER.info("Reading patient " + patient.patientId());
        String gender = null;
        String ethnicity = null;
        String tumorLocation = null;
        LocalDate registrationDate = null;
        String hospital1 = null;
        String hospital2 = null;
        String birthYear1 = null;
        String birthYear2 = null;
        LocalDate birthYear3 = null;
        LocalDate deathDate = null;
        final String impliedHospital = getHospital(patient, hospitals);

        for (final EcrfStudyEvent studyEvent : patient.studyEventsPerOID(STUDY_BASELINE)) {
            for (final EcrfForm demographyForm : studyEvent.nonEmptyFormsPerOID(FORM_DEMOGRAPHY, true)) {
                for (final EcrfItemGroup demographyItemGroup : demographyForm.nonEmptyItemGroupsPerOID(
                        ITEMGROUP_DEMOGRAPHY, true)) {
                    gender = demographyItemGroup.readItemString(FIELD_SEX, 0, true);
                    ethnicity = demographyItemGroup.readItemString(FIELD_ETHNICITY, 0, true);
                }
            }

            for (final EcrfForm carcinomaForm : studyEvent.nonEmptyFormsPerOID(FORM_CARCINOMA, true)) {
                for (final EcrfItemGroup carcinomaItemGroup : carcinomaForm.nonEmptyItemGroupsPerOID(
                        ITEMGROUP_CARCINOMA, true)) {
                    tumorLocation = carcinomaItemGroup.readItemString(FIELD_TUMOR_LOCATION, 0, true);
                }
            }

            for (final EcrfForm eligibilityForm : studyEvent.nonEmptyFormsPerOID(FORM_ELIGIBILITY, true)) {
                for (final EcrfItemGroup eligibilityItemGroup : eligibilityForm.nonEmptyItemGroupsPerOID(
                        ITEMGROUP_ELIGIBILITY, true)) {
                    registrationDate = eligibilityItemGroup.readItemDate(FIELD_REGISTRATION_DATE, 0, dateFormatter,
                            true);
                    hospital1 = eligibilityItemGroup.readItemString(FIELD_HOSPITAL1, 0, false);
                    birthYear2 = eligibilityItemGroup.readItemString(FIELD_BIRTH_YEAR2, 0, false);
                    birthYear3 = eligibilityItemGroup.readItemDate(FIELD_BIRTH_YEAR3, 0, dateFormatter, false);
                }
            }

            for (final EcrfForm selcritForm : studyEvent.nonEmptyFormsPerOID(FORM_SELCRIT, true)) {
                for (final EcrfItemGroup selcritItemGroup : selcritForm.nonEmptyItemGroupsPerOID(ITEMGROUP_SELCRIT,
                        true)) {
                    hospital2 = selcritItemGroup.readItemString(FIELD_HOSPITAL2, 0, false);
                    birthYear1 = selcritItemGroup.readItemString(FIELD_BIRTH_YEAR1, 0, false);
                }
            }
        }
        for (final EcrfStudyEvent endStudyEvent : patient.studyEventsPerOID(STUDY_ENDSTUDY)) {
            for (final EcrfForm deathFrom : endStudyEvent.nonEmptyFormsPerOID(FORM_DEATH, false)) {
                for (final EcrfItemGroup deathItemGroup : deathFrom.nonEmptyItemGroupsPerOID(ITEMGROUP_DEATH, false)) {
                    deathDate = deathItemGroup.readItemDate(FIELD_DEATH_DATE, 0, dateFormatter, true);
                }
            }
        }
        final Integer birthYear = determineBirthYear(birthYear1, birthYear2, birthYear3);
        if (birthYear == null) {
            LOGGER.warn(patient.patientId() + ": empty birth year fields");
        }
        checkHospitalVsImplied(patient.patientId(), impliedHospital, hospital1, FIELD_HOSPITAL1);
        checkHospitalVsImplied(patient.patientId(), impliedHospital, hospital2, FIELD_HOSPITAL2);
        return new PatientInfo(patient.patientId(), registrationDate, gender, ethnicity, impliedHospital, birthYear,
                tumorLocation, deathDate);
    }

    @NotNull
    private static Map<Integer, String> getHospitals(@NotNull final CpctEcrfModel datamodel) {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        final Iterable<EcrfField> fields = datamodel.findFieldsById(
                Lists.newArrayList(DATAMODEL_HOSPITAL1, DATAMODEL_HOSPITAL2));
        fields.forEach(field -> hospitals.putAll(field.codeList()));
        return ImmutableMap.copyOf(hospitals);
    }

    @Nullable
    private static String getHospital(@NotNull final EcrfPatient patient,
            @NotNull final Map<Integer, String> hospitals) {
        final Integer hospitalCode = Integer.parseInt(patient.patientId().substring(6, 8));
        final String hospital = hospitals.get(hospitalCode);
        if (hospital == null) {
            LOGGER.warn(DATAMODEL_HOSPITAL1 + ", " + DATAMODEL_HOSPITAL2 + " contained no Hospital with code "
                    + hospitalCode);
        }
        return hospital;
    }

    private void checkHospitalVsImplied(@NotNull final String patientId, @Nullable final String impliedHospital,
            @Nullable final String hospital, @Nullable final String hospitalField) {
        if (impliedHospital != null && hospital != null && !hospital.equals(impliedHospital)) {
            LOGGER.warn(patientId + ": " + hospitalField + " value(" + hospital + ") does not match cpctId value("
                    + impliedHospital + ")");
        }
    }

    @Nullable
    private Integer determineBirthYear(@Nullable final String birthYear1, @Nullable final String birthYear2,
            @Nullable final LocalDate birthYear3) {
        if (birthYear1 != null) {
            return Integer.parseInt(birthYear1);
        }
        if (birthYear2 != null) {
            return Integer.parseInt(birthYear2);
        }
        if (birthYear3 != null) {
            return birthYear3.getYear();
        }
        return null;
    }
}
