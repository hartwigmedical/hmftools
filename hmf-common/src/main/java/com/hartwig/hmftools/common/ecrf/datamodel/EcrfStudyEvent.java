package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class EcrfStudyEvent {
    @NotNull
    private final Map<String, List<EcrfForm>> formsPerOID;

    public EcrfStudyEvent() {
        this.formsPerOID = Maps.newHashMap();
    }

    public void addForm(@NotNull final String formOid, @NotNull final EcrfForm form) {
        if (!formsPerOID.containsKey(formOid)) {
            formsPerOID.put(formOid, Lists.newArrayList());
        }
        formsPerOID.get(formOid).add(form);
    }

    @NotNull
    public Map<String, List<EcrfForm>> formsPerOID() {
        return formsPerOID;
    }

    @NotNull
    public List<EcrfForm> formsPerOID(@NotNull final String formOID) {
        if (formsPerOID.get(formOID) == null) {
            return Lists.newArrayList();
        }
        return formsPerOID.get(formOID);
    }
}
