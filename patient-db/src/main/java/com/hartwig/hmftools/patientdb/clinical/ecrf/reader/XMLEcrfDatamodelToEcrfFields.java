package com.hartwig.hmftools.patientdb.clinical.ecrf.reader;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.patientdb.clinical.ecrf.datamodel.ImmutableEcrfDatamodelField;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfDatamodelToEcrfFields {

    private XMLEcrfDatamodelToEcrfFields() {
    }

    @NotNull
    public static List<EcrfDatamodelField> convert(@NotNull final XMLEcrfDatamodel datamodel) {
        List<EcrfDatamodelField> fields = Lists.newArrayList();
        for (StudyEvent studyEvent : datamodel.studyEvents().values()) {
            for (String formOID : studyEvent.formOIDs()) {
                Form form = findByOID(datamodel.forms().values(), formOID);
                for (String itemGroupOID : form.itemGroupOIDs()) {
                    ItemGroup itemGroup = findByOID(datamodel.itemGroups().values(), itemGroupOID);
                    for (String itemOID : itemGroup.itemOIDs()) {
                        Item item = findByOID(datamodel.items().values(), itemOID);
                        Map<Integer, String> codeList = Maps.newHashMap();
                        String codeListOID = item.codeListOID();
                        if (codeListOID != null) {
                            codeList = findByOID(datamodel.codeLists().values(), codeListOID).values();
                        }
                        EcrfDatamodelField field = new ImmutableEcrfDatamodelField(studyEvent.oid(),
                                form.oid(),
                                itemGroup.oid(),
                                item.oid(),
                                item.name(),
                                codeList);
                        fields.add(field);
                    }
                }
            }
        }
        return fields;
    }

    @NotNull
    private static <X extends OIDObject> X findByOID(@NotNull Collection<X> objects, @NotNull String OID) {
        for (X object : objects) {
            if (object.oid().equals(OID)) {
                return object;
            }
        }

        throw new IllegalStateException("Could not resolve " + OID);
    }
}
