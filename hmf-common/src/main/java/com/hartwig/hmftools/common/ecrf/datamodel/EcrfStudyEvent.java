package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class EcrfStudyEvent {
    private static final Logger LOGGER = LogManager.getLogger(EcrfForm.class);

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<String, List<EcrfForm>> formsPerOID;

    public EcrfStudyEvent(@NotNull final String patientId) {
        this.patientId = patientId;
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
    public List<EcrfForm> nonEmptyFormsPerOID(@NotNull final String formOID) {
        return nonEmptyFormsPerOID(formOID, false);
    }

    @NotNull
    public List<EcrfForm> nonEmptyFormsPerOID(@NotNull final String formOID, boolean verbose) {
        final List<EcrfForm> nonEmptyForms = Lists.newArrayList();
        if (formsPerOID.get(formOID) == null) {
            return Lists.newArrayList();
        } else {
            for (final EcrfForm form : formsPerOID.get(formOID)) {
                if (form.isEmpty() && verbose) {
                    LOGGER.warn(patientId + ": empty form: " + formOID);
                }
                if (!form.isEmpty()) {
                    nonEmptyForms.add(form);
                }
            }
        }
        return nonEmptyForms;
    }
}
