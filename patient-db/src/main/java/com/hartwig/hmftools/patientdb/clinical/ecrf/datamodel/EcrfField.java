package com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel;

import java.util.List;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public interface EcrfField {

    String OID_SEPARATOR = ".";

    List<String> IRRELEVANT_ITEM_GROUP_OIDS =
            Lists.newArrayList("AUDITTRAILENTRIES", "METADATA", "AUDITDATA", "QUERY", "VALIDATIONENTRIES");
    List<String> IRRELEVANT_ITEM_OIDS = Lists.newArrayList("AUDIT_ENTRIES", "FORMSTATUS");
    String IRRELEVANT_ITEM_OID_PATTERN = "GROUP";

    @NotNull
    String studyEventOID();

    @NotNull
    String formOID();

    @NotNull
    String itemGroupOID();

    @NotNull
    String itemOID();

    @NotNull
    default String name() {
        return EcrfFieldFunctions.name(studyEventOID(), formOID(), itemGroupOID(), itemOID());
    }

    default boolean isRelevant() {
        String[] parts = name().split("\\" + OID_SEPARATOR);
        String itemGroup = parts[2];
        String item = parts[3];
        return !(IRRELEVANT_ITEM_GROUP_OIDS.contains(itemGroup) || IRRELEVANT_ITEM_OIDS.contains(item) || item.contains(
                IRRELEVANT_ITEM_OID_PATTERN));
    }
}
