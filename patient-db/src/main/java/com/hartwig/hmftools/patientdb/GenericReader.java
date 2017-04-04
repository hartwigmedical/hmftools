package com.hartwig.hmftools.patientdb;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

class GenericReader {
    private static final Logger LOGGER = LogManager.getLogger(GenericReader.class);

    @Nullable
    static String getField(@NotNull final EcrfPatient patient, @NotNull final String fieldName) {
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
    static List<String> getFieldValues(@NotNull final EcrfPatient patient, @NotNull final String fieldName) {
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
        } else {
            for (int index = 0; index < values.size(); index++) {
                if (values.get(index) == null) {
                    LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains null at index " + index
                            + ".");
                } else if (values.get(index).replaceAll("\\s", "").length() == 0) {
                    LOGGER.warn(
                            fieldName + " for patient " + patient.patientId() + " contains only whitespaces at index "
                                    + index + ".");
                }
            }
        }
        return values;
    }

    @Nullable
    static List<String> getFieldValuesWithOthers(@NotNull final EcrfPatient patient, @NotNull final String fieldName,
            @NotNull final String othersFieldName) {
        final List<String> values = getFieldValues(patient, fieldName);
        if (hasOthers(values)) {
            final List<String> otherValues = getFieldValues(patient, othersFieldName);
            return replaceOthers(values, otherValues);
        } else {
            return values;
        }
    }

    private static boolean hasOthers(@Nullable final List<String> values) {
        if (values != null && values.size() > 0) {
            for (final String value : values) {
                if (value.toLowerCase().contains("other")) {
                    return true;
                }
            }
        }
        return false;
    }

    @Nullable
    private static List<String> replaceOthers(@NotNull final List<String> values,
            @Nullable final List<String> otherValues) {
        final List<String> result = Lists.newArrayList();
        for (int index = 0; index < values.size(); index++) {
            if (values.get(index).replaceAll("\\s", "").toLowerCase().startsWith("other")) {
                final String otherValue = Utils.getElemAtIndex(otherValues, index);
                result.add(otherValue);
            } else {
                result.add(values.get(index));
            }
        }
        return result;
    }

}
