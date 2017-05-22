package com.hartwig.hmftools.common.ecrf.datamodel;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfItemGroup {
    private static final Logger LOGGER = LogManager.getLogger(EcrfItemGroup.class);

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<String, List<String>> itemsPerOID;

    public EcrfItemGroup(@NotNull final String patientId) {
        this.patientId = patientId;
        this.itemsPerOID = Maps.newHashMap();
    }

    public void addItem(@NotNull final String itemOid, @NotNull final String itemValue) {
        if (!itemsPerOID.containsKey(itemOid)) {
            itemsPerOID.put(itemOid, Lists.newArrayList());
        }
        itemsPerOID.get(itemOid).add(itemValue);
    }

    @NotNull
    public Map<String, List<String>> itemsPerOID() {
        return itemsPerOID;
    }

    @NotNull
    public List<String> itemsPerOID(@NotNull final String itemOID) {
        if (itemsPerOID().get(itemOID) == null) {
            return Lists.newArrayList();
        }
        return itemsPerOID.get(itemOID);
    }

    public boolean isEmpty() {
        return itemsPerOID.values().stream().filter(
                listOfValues -> listOfValues.stream().filter(value -> value != null && !value.trim().isEmpty()).count()
                        > 0).count() == 0;
    }

    @Nullable
    private String readItemString(@NotNull final String itemOID, int index) {
        if (index < itemsPerOID(itemOID).size()) {
            final String ecrfValue = itemsPerOID(itemOID).get(index);
            if (ecrfValue != null && ecrfValue.replaceAll("\\s", "").length() == 0) {
                return null;
            }
            return ecrfValue;
        }
        return null;
    }

    @Nullable
    private LocalDate readItemDate(@NotNull final String itemOID, int index,
            @NotNull final DateTimeFormatter dateFormatter) {
        if (index < itemsPerOID(itemOID).size()) {
            final String ecrfValue = itemsPerOID(itemOID).get(index);
            if (ecrfValue == null) {
                return null;
            }
            try {
                return LocalDate.parse(ecrfValue, dateFormatter);
            } catch (DateTimeParseException e) {
                return null;
            }
        }
        return null;
    }

    @Nullable
    public LocalDate readItemDate(@NotNull final String itemOID, int index,
            @NotNull final DateTimeFormatter dateFormatter, boolean verbose) {
        final LocalDate itemDate = readItemDate(itemOID, index, dateFormatter);
        if (itemDate == null && verbose) {
            LOGGER.warn(patientId + ": empty field: " + itemOID);
        }
        return itemDate;
    }

    @Nullable
    public String readItemString(@NotNull final String itemOID, int index, boolean verbose) {
        final String itemString = readItemString(itemOID, index);
        if (itemString == null && verbose) {
            LOGGER.warn(patientId + ": empty field: " + itemOID);
        }
        return itemString;
    }
}
