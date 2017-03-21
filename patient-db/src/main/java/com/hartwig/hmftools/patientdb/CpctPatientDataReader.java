package com.hartwig.hmftools.patientdb;

import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.List;
import java.util.Map;

import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.CpctEcrfModel;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CpctPatientDataReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);
    private static final DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");
    private final CpctEcrfModel dataModel;
    private final Map<Integer, String> hospitals;

    public CpctPatientDataReader(@NotNull CpctEcrfModel model) {
        this.dataModel = model;
        this.hospitals = getHospitals(model);
    }

    @NotNull
    public PatientData read(@NotNull EcrfPatient patient) {
        final String sex = GenericReader.getField(patient, "BASELINE.DEMOGRAPHY.DEMOGRAPHY.SEX");
        final String ethnicity = GenericReader.getField(patient, "BASELINE.DEMOGRAPHY.DEMOGRAPHY.ETHNIC");
        final Integer birthYear = getBirthYear(patient, dateFormat);
        final String hospital = getHospital(patient, hospitals);
        final PatientData patientData = new PatientData(patient.patientId(), null, sex, birthYear, hospital,
                ethnicity);
        return patientData;
    }

    @NotNull
    private static Map<Integer, String> getHospitals(@NotNull CpctEcrfModel datamodel) {
        final Map<Integer, String> hospitals = Maps.newHashMap();
        final Iterable<EcrfField> fields = datamodel.findFieldsById(
                Lists.newArrayList("BASELINE.ELIGIBILITY.ELIGIBILITY.HOSPITAL", "BASELINE.SELCRIT.SELCRIT.NHOSPITAL"));
        fields.forEach(field -> hospitals.putAll(field.codeList()));
        return ImmutableMap.copyOf(hospitals);
    }

    @Nullable
    private static Integer getBirthYear(@NotNull EcrfPatient patient, @NotNull DateFormat dateFormat) {
        final List<String> nbirthYearValues = patient.fieldValuesByName("BASELINE.SELCRIT.SELCRIT.NBIRTHYEAR");
        if (nbirthYearValues != null && nbirthYearValues.size() != 0) {
            final String nbirthYear = nbirthYearValues.get(0);
            if (nbirthYear != null && !nbirthYear.equals("")) {
                return Integer.parseInt(nbirthYear);
            }
        }
        final List<String> birthYearValues = patient.fieldValuesByName("BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHYEAR");
        if (birthYearValues != null && birthYearValues.size() != 0) {
            final String birthYear = birthYearValues.get(0);
            if (birthYear != null && !birthYear.equals("")) {
                return Integer.parseInt(birthYear);
            }
        }
        final List<String> birthdtcesValues = patient.fieldValuesByName("BASELINE.ELIGIBILITY.ELIGIBILITY.BIRTHDTCES");
        if (birthdtcesValues != null && birthdtcesValues.size() != 0) {
            final String birthdtces = birthdtcesValues.get(0);
            if (birthdtces != null && !birthdtces.equals("")) {
                try {
                    final Date date = dateFormat.parse(birthdtces);
                    final Calendar calendar = Calendar.getInstance();
                    calendar.setTime(date);
                    return calendar.get(Calendar.YEAR);
                } catch (java.text.ParseException e) {
                    LOGGER.info("BIRTHDTCES field did not contain valid date for patient " + patient.patientId());
                    return null;
                }
            }
        }
        LOGGER.info("No completed birth year fields found for patient " + patient.patientId());
        return null;
    }

    @NotNull
    private static String getHospital(@NotNull EcrfPatient patient, @NotNull Map<Integer, String> hospitals) {
        final Integer hospitalCode = Integer.parseInt(patient.patientId().substring(6, 8));
        return hospitals.get(hospitalCode);
    }
}
