package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfDatamodelToEcrfFields {

    private XMLEcrfDatamodelToEcrfFields() {
    }

    @NotNull
    public static List<EcrfField> convert(@NotNull XMLEcrfDatamodel datamodel) {
        List<EcrfField> fields = Lists.newArrayList();
        for (StudyEvent studyEvent : datamodel.studyEvents()) {
            for (String formOID : studyEvent.formOIDs()) {
                Form form = findByOID(datamodel.forms(), formOID);
                for (String itemGroupOID : form.itemGroupOIDs()) {
                    ItemGroup itemGroup = findByOID(datamodel.itemGroups(), itemGroupOID);
                    for (String itemOID : itemGroup.itemOIDs()) {
                        Item item = findByOID(datamodel.items(), itemOID);
                        Map<Integer, String> codeList = Maps.newHashMap();
                        String codeListOID = item.codeListOID();
                        if (codeListOID != null) {
                            CodeList codeListObj = findByOID(datamodel.codeLists(), codeListOID);
                            codeList = codeListObj.values();
                        }
                        //                        String name = OIDFunctions.toEcrfFieldName(item.OID());
                        fields.add(
                                new EcrfField(studyEvent.OID(), form.OID(), itemGroup.OID(), item.OID(), item.name(),
                                        codeList));
                    }
                }
            }
        }
        return fields;
    }

    @NotNull
    private static <X extends OIDObject> X findByOID(@NotNull List<X> objects, @NotNull String OID) {
        for (X object : objects) {
            if (object.OID().equals(OID)) {
                return object;
            }
        }

        throw new IllegalStateException("Could not resolve " + OID + " as a " + objects.get(0).getClass().getName());
    }
}
