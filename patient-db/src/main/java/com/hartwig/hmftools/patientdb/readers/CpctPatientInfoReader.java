package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.Utils;
import com.hartwig.hmftools.patientdb.data.PatientInfo;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CpctPatientInfoReader {
    private static final Logger LOGGER = LogManager.getLogger(CpctPatientInfoReader.class);

    private static final String FIELD_SEX = "BASELINE.DEMOGRAPHY.DEMOGRAPHY.SEX";
    private static final String FIELD_ETHNICITY = "BASELINE.DEMOGRAPHY.DEMOGRAPHY.ETHNIC";
    private static final String FIELD_HOSPITALS1 = "BASELINE.ELIGIBILITY.ELIGIBILITY.HOSPITAL";
    private static final String FIELD_HOSPITALS2 = "BASELINE.SELCRIT.SELCRIT.NHOSPITAL";
    private static final String FIELD_BIRTHYEAR1 = "BASELINE.SELCRIT.SELCRIT.NBIRTHYEAR";
    private static final String FIELD_BIRTHYEAR2 = "BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHYEAR";
    private static final String FIELD_BIRTHYEAR3 = "BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHDTCES";
    private static final String FIELD_DEATHDATE = "ENDSTUDY.DEATH.DEATH.DDEATHDTC";

    private static final DateTimeFormatter dateFormatter = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final Map<Integer, String> hospitals;

    public CpctPatientInfoReader(@NotNull final CpctEcrfModel model) {
        this.hospitals = getHospitals(model);
    }

    @NotNull
    public PatientInfo read(@NotNull final EcrfPatient patient) {
        final String sex = GenericReader.getField(patient, FIELD_SEX);
        final String ethnicity = GenericReader.getField(patient, FIELD_ETHNICITY);
        final Integer birthYear = getBirthYear(patient);
        final String hospital = getHospital(patient, hospitals);
        final String deathDateString = GenericReader.getField(patient, FIELD_DEATHDATE);
        final LocalDate deathDate = Utils.getDate(deathDateString, dateFormatter);
        return new PatientInfo(patient.patientId(), null, sex, birthYear, hospital, ethnicity, deathDate);
    }

    @NotNull
    private static Map<Integer, String> getHospitals(@NotNull final CpctEcrfModel datamodel) {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        final Iterable<EcrfField> fields = datamodel.findFieldsById(
                Lists.newArrayList(FIELD_HOSPITALS1, FIELD_HOSPITALS2));
        fields.forEach(field -> hospitals.putAll(field.codeList()));
        return ImmutableMap.copyOf(hospitals);
    }

    @Nullable
    private static Integer getBirthYear(@NotNull final EcrfPatient patient) {
        final List<String> nbirthYearValues = patient.fieldValuesByName(FIELD_BIRTHYEAR1);
        if (nbirthYearValues != null && nbirthYearValues.size() != 0) {
            final String nbirthYear = nbirthYearValues.get(0);
            if (nbirthYear != null && !nbirthYear.equals("")) {
                return Integer.parseInt(nbirthYear);
            }
        }
        final List<String> birthYearValues = patient.fieldValuesByName(FIELD_BIRTHYEAR2);
        if (birthYearValues != null && birthYearValues.size() != 0) {
            final String birthYear = birthYearValues.get(0);
            if (birthYear != null && !birthYear.equals("")) {
                return Integer.parseInt(birthYear);
            }
        }
        final List<String> birthdtcesValues = patient.fieldValuesByName(FIELD_BIRTHYEAR3);
        if (birthdtcesValues != null && birthdtcesValues.size() != 0) {
            final String birthdtces = birthdtcesValues.get(0);
            if (birthdtces != null && !birthdtces.equals("")) {
                try {
                    final LocalDate date = LocalDate.parse(birthdtces, dateFormatter);
                    return date.getYear();
                } catch (DateTimeParseException e) {
                    LOGGER.warn(FIELD_BIRTHYEAR3 + " did not contain valid date for patient " + patient.patientId());
                    return null;
                }
            }
        }
        LOGGER.warn("No completed birth year fields found for patient " + patient.patientId());
        return null;
    }

    @Nullable
    private static String getHospital(@NotNull final EcrfPatient patient,
            @NotNull final Map<Integer, String> hospitals) {
        final Integer hospitalCode = Integer.parseInt(patient.patientId().substring(6, 8));
        final String hospital = hospitals.get(hospitalCode);
        if (hospital == null) {
            LOGGER.warn(
                    FIELD_HOSPITALS1 + ", " + FIELD_HOSPITALS2 + " contained no Hospital with code " + hospitalCode);
        }
        return hospital;
    }
}
