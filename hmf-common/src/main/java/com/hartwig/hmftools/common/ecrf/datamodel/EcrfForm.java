package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class EcrfForm {
    @NotNull
    private final Map<String, List<EcrfItemGroup>> itemGroupsPerOID;

    public EcrfForm() {
        this.itemGroupsPerOID = Maps.newHashMap();
    }

    public void addItemGroup(@NotNull final String itemGroupOid, @NotNull final EcrfItemGroup itemGroup) {
        if (!itemGroupsPerOID.containsKey(itemGroupOid)) {
            itemGroupsPerOID.put(itemGroupOid, Lists.newArrayList());
        }
        itemGroupsPerOID.get(itemGroupOid).add(itemGroup);
    }

    @NotNull
    public Map<String, List<EcrfItemGroup>> itemGroupsPerOID() {
        return itemGroupsPerOID;
    }

    @NotNull
    public List<EcrfItemGroup> itemGroupsPerOID(@NotNull final String itemGroupOID) {
        if (itemGroupsPerOID.get(itemGroupOID) == null) {
            return Lists.newArrayList();
        }
        return itemGroupsPerOID.get(itemGroupOID);
    }

    public boolean isEmpty() {
        return itemGroupsPerOID.values().stream().filter(
                itemGroups -> itemGroups.stream().filter(group -> !group.isEmpty()).count() > 0).count() == 0;
    }

}
