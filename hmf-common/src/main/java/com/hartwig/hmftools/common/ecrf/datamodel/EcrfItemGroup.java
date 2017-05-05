package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class EcrfItemGroup {
    @NotNull
    private final Map<String, List<String>> itemsPerOID;

    public EcrfItemGroup() {
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
}
