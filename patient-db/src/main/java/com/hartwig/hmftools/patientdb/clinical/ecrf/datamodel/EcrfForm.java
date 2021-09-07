package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.ecrf.formstatus.FormStatus;

import org.jetbrains.annotations.NotNull;

public class EcrfForm {

    @NotNull
    private final FormStatus formStatus;
    @NotNull
    private final Map<String, List<EcrfItemGroup>> itemGroupsPerOID;

    public EcrfForm(@NotNull final FormStatus formStatus) {
        this.formStatus = formStatus;
        this.itemGroupsPerOID = Maps.newHashMap();
    }

    public void addItemGroup(@NotNull String itemGroupOid, @NotNull EcrfItemGroup itemGroup) {
        if (!itemGroupsPerOID.containsKey(itemGroupOid)) {
            itemGroupsPerOID.put(itemGroupOid, Lists.newArrayList());
        }
        itemGroupsPerOID.get(itemGroupOid).add(itemGroup);
    }

    @NotNull
    public Map<String, List<EcrfItemGroup>> itemGroupsPerOID() {
        return itemGroupsPerOID;
    }

    public boolean isEmpty() {
        return itemGroupsPerOID.values().stream().noneMatch(itemGroups -> itemGroups.stream().anyMatch(group -> !group.isEmpty()));
    }

    @NotNull
    public List<EcrfItemGroup> nonEmptyItemGroupsPerOID(@NotNull String itemGroupOID) {
        List<EcrfItemGroup> nonEmptyItemGroups = Lists.newArrayList();
        if (itemGroupsPerOID.get(itemGroupOID) == null) {
            return Lists.newArrayList();
        } else {
            for (EcrfItemGroup itemGroup : itemGroupsPerOID.get(itemGroupOID)) {
                if (!itemGroup.isEmpty()) {
                    nonEmptyItemGroups.add(itemGroup);
                }
            }
        }
        return nonEmptyItemGroups;
    }

    public boolean locked() {
        return formStatus.locked();
    }

    @NotNull
    public FormStatus status() {
        return formStatus;
    }
}
