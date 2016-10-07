package com.hartwig.hmftools.ecrfanalyser.datamodel;

import java.util.Map;

import org.jetbrains.annotations.NotNull;

public class EcrfField implements Comparable<EcrfField> {

    private static final String BIRTH_DATE_IDENTIFIER = "BIRTHDTC";

    @NotNull
    private final String studyEventOID;
    @NotNull
    private final String formOID;
    @NotNull
    private final String itemGroupOID;
    @NotNull
    private final String itemOID;
    @NotNull
    private final String description;
    @NotNull
    private final Map<Integer, String> codeList;

    public EcrfField(@NotNull final String studyEventOID, @NotNull final String formOID,
            @NotNull final String itemGroupOID, @NotNull final String itemOID, @NotNull final String description,
            @NotNull final Map<Integer, String> codeList) {
        this.studyEventOID = studyEventOID;
        this.formOID = formOID;
        this.itemGroupOID = itemGroupOID;
        this.itemOID = itemOID;
        this.description = description;
        this.codeList = codeList;
    }

    @NotNull
    public String name() {
        return itemOID;
    }

    @NotNull
    public String studyEventOID() {
        return studyEventOID;
    }

    @NotNull
    public String formOID() {
        return formOID;
    }

    @NotNull
    public String itemGroupOID() {
        return itemGroupOID;
    }

    @NotNull
    public String itemOID() {
        return itemOID;
    }

    @NotNull
    public String description() {
        return description;
    }

    @NotNull
    public Map<Integer, String> codeList() {
        return codeList;
    }

    @NotNull
    public String resolveValue(@NotNull String ecrfValue) throws EcrfResolveException {
        String value;
        if (codeList.size() > 0 && ecrfValue.length() > 0) {
            if (isInteger(ecrfValue)) {
                value = codeList.get(Integer.valueOf(ecrfValue));
                if (value == null) {
                    throw new EcrfResolveException(
                            "Could not find value in dropdown for " + itemOID + ": " + ecrfValue);
                }
            } else {
                throw new EcrfResolveException(
                        "Could not convert value from list to integer for " + itemOID + ": " + ecrfValue);
            }
        } else {
            value = ecrfValue;
        }

        if (itemOID.contains(BIRTH_DATE_IDENTIFIER) && value.length() > 0) {
            if (value.length() < 4) {
                throw new EcrfResolveException("Could not convert " + name() + ": " + value);
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

    public int compareTo(@NotNull EcrfField other) {
        return name().compareTo(other.name());
    }
}
