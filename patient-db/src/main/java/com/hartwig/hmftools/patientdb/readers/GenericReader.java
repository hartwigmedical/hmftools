package com.hartwig.hmftools.patientdb.readers;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfPatient;
import com.hartwig.hmftools.patientdb.Utils;

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
        final List<String> result = Lists.newArrayList();
        final List<String> values = patient.fieldValuesByName(fieldName);
        if (values == null) {
            LOGGER.warn(fieldName + " not found for patient " + patient.patientId() + ".");
            return null;
        } else if (values.size() == 0) {
            LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains no values.");
        } else {
            for (int index = 0; index < values.size(); index++) {
                final String fieldValue = values.get(index);
                if (fieldValue == null) {
                    LOGGER.warn(fieldName + " for patient " + patient.patientId() + " contains null at index " + index
                            + ".");
                    result.add(null);
                } else if (fieldValue.replaceAll("\\s", "").length() == 0) {
                    LOGGER.warn(
                            fieldName + " for patient " + patient.patientId() + " contains only whitespaces at index "
                                    + index + ".");
                    result.add(null);
                } else {
                    result.add(fieldValue);
                }
            }
        }
        return result;
    }

    @Nullable
    static List<LocalDate> getDateFieldValues(@NotNull final EcrfPatient patient, @NotNull final String fieldName,
            @NotNull final DateTimeFormatter dateFormatter) {
        final List<LocalDate> fieldDates = Lists.newArrayList();
        final List<String> fieldValues = getFieldValues(patient, fieldName);
        for (int index = 0; index < fieldValues.size(); index++) {
            final String fieldValue = fieldValues.get(index);
            final LocalDate fieldDate = Utils.getDate(fieldValue, dateFormatter);
            if (fieldDate == null) {
                LOGGER.warn(fieldName + " did not contain valid date at index " + index + " for patient  "
                        + patient.patientId());
            }
            fieldDates.add(fieldDate);
        }
        return fieldDates;
    }

    /**
     * Reads the values for fieldName and replaces values that start with 'other' with values from othersFieldName.
     * Meant to be used for fields like ethnicity, where the last possible value is 'other or unknown', in which case
     * a different (free text) field is used to fill in ethnicity data.
     *
     * @param patient         patient to read fields from.
     * @param fieldName       name of field to read data form.
     * @param othersFieldName name of field to read data in case the value in fieldName is something like "other, please specify"
     * @return list of values for fieldName with 'other' values replaced from othersFieldName
     */
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
                if (value != null && value.replaceAll("\\s", "").toLowerCase().startsWith("other")) {
                    return true;
                }
            }
        }
        return false;
    }

    /**
     * Iterates over the values list and replaces any value that starts with 'other' with the item found
     * at the same index in the otherValues list.
     *
     * @param values      list of values to check
     * @param otherValues list of replacement values
     * @return list of items from the values list, with items starting with 'other' replaced
     */
    @Nullable
    private static List<String> replaceOthers(@NotNull final List<String> values,
            @Nullable final List<String> otherValues) {
        final List<String> result = Lists.newArrayList();
        for (int index = 0; index < values.size(); index++) {
            final String value = values.get(index);
            if (value != null && value.replaceAll("\\s", "").toLowerCase().startsWith("other")) {
                final String otherValue = Utils.getElemAtIndex(otherValues, index);
                result.add(otherValue);
            } else {
                result.add(value);
            }
        }
        return result;
    }
}
