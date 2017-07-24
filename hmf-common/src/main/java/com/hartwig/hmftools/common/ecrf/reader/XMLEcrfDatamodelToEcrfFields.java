package com.hartwig.hmftools.common.ecrf.reader;

import java.util.Collection;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.ecrf.datamodel.EcrfDatamodelField;
import com.hartwig.hmftools.common.ecrf.datamodel.ImmutableEcrfDatamodelField;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfDatamodelToEcrfFields {

    private XMLEcrfDatamodelToEcrfFields() {
    }

    @NotNull
    public static List<EcrfDatamodelField> convert(@NotNull final XMLEcrfDatamodel datamodel) {
        final List<EcrfDatamodelField> fields = Lists.newArrayList();
        for (final StudyEvent studyEvent : datamodel.studyEvents().values()) {
            for (final String formOID : studyEvent.formOIDs()) {
                final Form form = findByOID(datamodel.forms().values(), formOID);
                for (final String itemGroupOID : form.itemGroupOIDs()) {
                    final ItemGroup itemGroup = findByOID(datamodel.itemGroups().values(), itemGroupOID);
                    for (final String itemOID : itemGroup.itemOIDs()) {
                        final Item item = findByOID(datamodel.items().values(), itemOID);
                        Map<Integer, String> codeList = Maps.newHashMap();
                        String codeListOID = item.codeListOID();
                        if (codeListOID != null) {
                            final CodeList codeListObj = findByOID(datamodel.codeLists().values(), codeListOID);
                            codeList = codeListObj.values();
                        }
                        final EcrfDatamodelField field =
                                new ImmutableEcrfDatamodelField(studyEvent.OID(), form.OID(), itemGroup.OID(), item.OID(), item.name(),
                                        codeList);
                        fields.add(field);
                    }
                }
            }
        }
        return fields;
    }

    @NotNull
    private static <X extends OIDObject> X findByOID(@NotNull final Collection<X> objects, @NotNull final String OID) {
        for (final X object : objects) {
            if (object.OID().equals(OID)) {
                return object;
            }
        }

        throw new IllegalStateException("Could not resolve " + OID);
    }
}
