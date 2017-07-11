package com.hartwig.hmftools.common.ecrf.datamodel;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class EcrfPatient {

    @NotNull
    private final String patientId;
    @NotNull
    private final Map<EcrfField, List<String>> fields;
    @NotNull
    private final Map<String, List<EcrfStudyEvent>> studyEventsPerOID;

    public EcrfPatient(@NotNull final String patientId, @NotNull final Map<String, List<EcrfStudyEvent>> studyEventsPerOID) {
        this.patientId = patientId;
        this.studyEventsPerOID = studyEventsPerOID;
        this.fields = createFields(studyEventsPerOID);
    }

    @NotNull
    public String patientId() {
        return patientId;
    }

    @NotNull
    public Set<EcrfField> fields() {
        return fields.keySet();
    }

    @Nullable
    public List<String> fieldValuesByEcrfField(@NotNull final EcrfField field) {
        return fields.get(field);
    }

    @Nullable
    public List<String> fieldValuesByName(@NotNull final String fieldName) {
        for (final EcrfField field : fields.keySet()) {
            if (field.name().equals(fieldName)) {
                return fieldValuesByEcrfField(field);
            }
        }
        return null;
    }

    @NotNull
    public Map<String, List<EcrfStudyEvent>> studyEventsPerOID() {
        return studyEventsPerOID;
    }

    @NotNull
    public List<EcrfStudyEvent> studyEventsPerOID(@NotNull final String studyEventOID) {
        if (studyEventsPerOID().get(studyEventOID) == null) {
            return Lists.newArrayList();
        }
        return studyEventsPerOID.get(studyEventOID);
    }

    @NotNull
    private static Map<EcrfField, List<String>> createFields(@NotNull final Map<String, List<EcrfStudyEvent>> studyEventsPerOID) {
        final Map<EcrfField, List<String>> fields = Maps.newHashMap();
        for (final String studyEventOID : studyEventsPerOID.keySet()) {
            for (final EcrfStudyEvent studyEvent : studyEventsPerOID.get(studyEventOID)) {
                for (final String formOID : studyEvent.formsPerOID().keySet()) {
                    for (final EcrfForm form : studyEvent.formsPerOID().get(formOID)) {
                        for (final String itemGroupOID : form.itemGroupsPerOID().keySet()) {
                            for (final EcrfItemGroup itemGroup : form.itemGroupsPerOID().get(itemGroupOID)) {
                                for (final String itemOID : itemGroup.itemsPerOID().keySet()) {
                                    final EcrfField field =
                                            new ImmutableEcrfField(studyEventOID, formOID, itemGroupOID, itemOID, "", Maps.newHashMap());
                                    fields.put(field, itemGroup.itemsPerOID().get(itemOID));
                                }
                            }
                        }
                    }
                }
            }
        }
        return fields;
    }
}
