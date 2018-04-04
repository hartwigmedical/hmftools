package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;

public class EcrfStudyEvent {

    @NotNull
    private final Map<String, List<EcrfForm>> formsPerOID = Maps.newHashMap();

    public EcrfStudyEvent() {
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
    public List<EcrfForm> nonEmptyFormsPerOID(@NotNull final String formOID) {
        final List<EcrfForm> nonEmptyForms = Lists.newArrayList();
        if (formsPerOID.get(formOID) == null) {
            return Lists.newArrayList();
        } else {
            for (final EcrfForm form : formsPerOID.get(formOID)) {
                if (!form.isEmpty()) {
                    nonEmptyForms.add(form);
                }
            }
        }
        return nonEmptyForms;
    }
}
