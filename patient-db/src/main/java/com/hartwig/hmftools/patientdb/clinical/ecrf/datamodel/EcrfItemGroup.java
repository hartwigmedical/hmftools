package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import java.time.LocalDate;
import java.time.format.DateTimeFormatter;
import java.time.format.DateTimeParseException;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfItemGroup {

    private static final DateTimeFormatter DATE_FORMATTER = DateTimeFormatter.ofPattern("yyyy-MM-dd");

    @NotNull
    private final Map<String, List<String>> itemsPerOID = Maps.newHashMap();

    public EcrfItemGroup() {
    }

    public void addItem(@NotNull String itemOid, @NotNull String itemValue) {
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
    private List<String> itemsPerOID(@NotNull String itemOID) {
        if (itemsPerOID().get(itemOID) == null) {
            return Lists.newArrayList();
        }
        return itemsPerOID.get(itemOID);
    }

    public boolean isEmpty() {
        return itemsPerOID.values()
                .stream()
                .noneMatch(listOfValues -> listOfValues.stream().anyMatch(value -> value != null && !value.trim().isEmpty()));
    }

    @Nullable
    public LocalDate readItemDate(@NotNull String itemOID) {
        String ecrfValue = readItemString(itemOID);

        if (ecrfValue == null) {
            return null;
        }

        try {
            return LocalDate.parse(ecrfValue.trim(), DATE_FORMATTER);
        } catch (DateTimeParseException e) {
            return null;
        }
    }

    @Nullable
    public String readItemString(@NotNull String itemOID) {
        // In theory we could have multiple values for one itemOID in one item group but in practice there is always a 1:1 relation.
        if (!itemsPerOID(itemOID).isEmpty()) {
            String ecrfValue = itemsPerOID(itemOID).get(0);
            if (ecrfValue != null) {
                if (ecrfValue.replaceAll("\\s", "").length() == 0) {
                    return null;
                } else {
                    // Remove whitespace + non-breakable spaces
                    return ecrfValue.trim().replaceAll("\\u00A0", "");
                }
            }
        }
        return null;
    }
}
