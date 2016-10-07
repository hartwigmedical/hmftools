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
    public static List<EcrfField> convert(@NotNull XMLEcrfDatamodel container) {
        List<EcrfField> fields = Lists.newArrayListWithCapacity(container.items().size());
        for (Item item : container.items()) {
            Map<Integer, String> values = Maps.newHashMap();
            if (item.codeListOID() != null) {
                values = findValuesForCodeList(container.codeLists(), item.codeListOID());
            }
            String name = OIDFunctions.toEcrfFieldName(item.OID());
            fields.add(new EcrfField(name, item.name(), values));
        }
        return fields;
    }

    @NotNull
    private static Map<Integer, String> findValuesForCodeList(final List<CodeList> codeLists,
            final String codeListOID) {
        for (CodeList codeList : codeLists) {
            if (codeList.OID().equals(codeListOID)) {
                return codeList.values();
            }
        }

        throw new IllegalStateException("Could not find code list : " + codeListOID);
    }

}
