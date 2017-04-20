package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfField;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfFieldFunctions;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfDatamodelToEcrfFields {

    private XMLEcrfDatamodelToEcrfFields() {
    }

    @NotNull
    public static List<EcrfField> convert(@NotNull final XMLEcrfDatamodel datamodel) {
        final List<EcrfField> fields = Lists.newArrayList();
        for (final StudyEvent studyEvent : datamodel.studyEvents()) {
            for (final String formOID : studyEvent.formOIDs()) {
                final Form form = findByOID(datamodel.forms(), formOID);
                for (final String itemGroupOID : form.itemGroupOIDs()) {
                    final ItemGroup itemGroup = findByOID(datamodel.itemGroups(), itemGroupOID);
                    for (final String itemOID : itemGroup.itemOIDs()) {
                        final Item item = findByOID(datamodel.items(), itemOID);
                        Map<Integer, String> codeList = Maps.newHashMap();
                        String codeListOID = item.codeListOID();
                        if (codeListOID != null) {
                            final CodeList codeListObj = findByOID(datamodel.codeLists(), codeListOID);
                            codeList = codeListObj.values();
                        }
                        final EcrfField field = new EcrfField(studyEvent.OID(), form.OID(), itemGroup.OID(),
                                item.OID(), item.name(), codeList);
                        if (EcrfFieldFunctions.isRelevant(field)) {
                            fields.add(field);
                        }
                    }
                }
            }
        }
        return fields;
    }

    @NotNull
    private static <X extends OIDObject> X findByOID(@NotNull final List<X> objects, @NotNull final String OID) {
        for (final X object : objects) {
            if (object.OID().equals(OID)) {
                return object;
            }
        }

        throw new IllegalStateException("Could not resolve " + OID + " as a " + objects.get(0).getClass().getName());
    }
}
