package com.hartwig.hmftools.patientdb;

import java.util.List;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class GenericReader {
    private static final Logger LOGGER = LogManager.getLogger(PatientDbRunner.class);

    @Nullable
    static String getField(@NotNull EcrfPatient patient, @NotNull String fieldName) {
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
            return null;
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
            return null;
        } else if (values.size() > 1) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains more than 1 value.");
        } else if (values.get(0).replaceAll("\\s", "").length() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains only whitespaces.");
            return null;
        }
        return values.get(0);
    }

    @Nullable
    static List<String> getFieldValues(@NotNull EcrfPatient patient, @NotNull String fieldName) {
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
            return null;
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
        }
        for (int index = 0; index < values.size(); index++) {
            if (values.get(index) == null) {
                LOGGER.warn(
                        fieldName + " for patient " + patient.patientId() + " contains null at index " + index + ".");
            } else if (values.get(index).replaceAll("\\s", "").length() == 0) {
                LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains only whitespaces at index "
                        + index + ".");
            }
        }
        return values;
    }
}
