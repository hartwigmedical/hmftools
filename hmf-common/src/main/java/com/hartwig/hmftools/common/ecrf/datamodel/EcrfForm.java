package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatus;
import com.hartwig.hmftools.common.ecrf.formstatus.FormStatusState;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfForm {
    private static final Logger LOGGER = LogManager.getLogger(EcrfForm.class);

    @NotNull
    private final String patientId;
    @NotNull
    private final FormStatus formStatus;

    @NotNull
    private final Map<String, List<EcrfItemGroup>> itemGroupsPerOID;

    public EcrfForm(@NotNull final String patientId, @NotNull final FormStatus formStatus) {
        this.patientId = patientId;
        this.formStatus = formStatus;
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

    public boolean isEmpty() {
        return itemGroupsPerOID.values().stream().noneMatch(itemGroups -> itemGroups.stream().anyMatch(group -> !group.isEmpty()));
    }

    @NotNull
    public List<EcrfItemGroup> nonEmptyItemGroupsPerOID(@NotNull final String itemGroupOID) {
        return nonEmptyItemGroupsPerOID(itemGroupOID, false);
    }

    @NotNull
    public List<EcrfItemGroup> nonEmptyItemGroupsPerOID(@NotNull final String itemGroupOID, boolean verbose) {
        final List<EcrfItemGroup> nonEmptyItemGroups = Lists.newArrayList();
        if (itemGroupsPerOID.get(itemGroupOID) == null) {
            return Lists.newArrayList();
        } else {
            for (final EcrfItemGroup itemGroup : itemGroupsPerOID.get(itemGroupOID)) {
                if (itemGroup.isEmpty() && verbose) {
                    LOGGER.warn(patientId + ": empty item group: " + itemGroupOID);
                }
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
    public FormStatusState status() {
        return formStatus.state();
    }
}
