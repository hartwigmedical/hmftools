package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class EcrfFieldFunctions {

    private static final String BIRTH_DATE_IDENTIFIER = "BIRTHDTC";
    private static final String OID_SEPARATOR = ".";

    private static final List<String> IRRELEVANT_ITEM_GROUP_OIDS = Lists.newArrayList("AUDITTRAILENTRIES", "METADATA",
            "AUDITDATA", "QUERY", "VALIDATIONENTRIES");
    private static final List<String> IRRELEVANT_ITEM_OIDS = Lists.newArrayList("AUDIT_ENTRIES", "FORMSTATUS");
    private static final String IRRELEVANT_ITEM_OID_PATTERN = "GROUP";

    private EcrfFieldFunctions() {
    }

    public static boolean isRelevant(@NotNull EcrfField field) {
        return isRelevant(field.name());
    }

    public static boolean isRelevant(@NotNull String name) {
        String parts[] = name.split("\\" + OID_SEPARATOR);
        String itemGroup = parts[2];
        String item = parts[3];
        return !(IRRELEVANT_ITEM_GROUP_OIDS.contains(itemGroup) || IRRELEVANT_ITEM_OIDS.contains(item)
                || item.contains(IRRELEVANT_ITEM_OID_PATTERN));
    }

    @NotNull
    public static String name(@NotNull String studyEventOID, @NotNull String formOID, @NotNull String itemGroupOID,
            @NotNull String itemOID) {
        String study = lastElement(studyEventOID);
        String form = lastElement(formOID);
        String itemGroup = lastElement(itemGroupOID);
        String item = lastElement(itemOID);

        return (study + OID_SEPARATOR + form + OID_SEPARATOR + itemGroup + OID_SEPARATOR + item).toUpperCase();
    }

    @NotNull
    private static String lastElement(@NotNull String OID) {
        String[] elements = OID.split("\\" + OID_SEPARATOR);
        return elements[elements.length - 1];
    }

    @NotNull
    public static String resolveValue(@NotNull EcrfField field, @NotNull String ecrfValue)
            throws EcrfResolveException {
        String value;
        if (field.codeList().size() > 0 && ecrfValue.length() > 0) {
            if (isInteger(ecrfValue)) {
                value = field.codeList().get(Integer.valueOf(ecrfValue));
                if (value == null) {
                    throw new EcrfResolveException(
                            "Could not find value in dropdown for " + field.name() + ": " + ecrfValue);
                }
            } else {
                throw new EcrfResolveException(
                        "Could not convert value from list to integer for " + field.name() + ": " + ecrfValue);
            }
        } else {
            value = ecrfValue;
        }

        if (field.itemOID().contains(BIRTH_DATE_IDENTIFIER) && value.length() > 0) {
            if (value.length() < 4) {
                throw new EcrfResolveException("Could not convert " + field.name() + ": " + value);
            } else {
                value = value.substring(0, 4) + "-01-01";
            }
        }

        return value;
    }

    private static boolean isInteger(@NotNull String integerString) {
        try {
            return Integer.valueOf(integerString) != null;
        } catch (NumberFormatException exception) {
            return false;
        }
    }
}
