package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;
import java.util.Map;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.stream.events.XMLEvent;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.ecrfanalyser.datamodel.EcrfField;

import org.jetbrains.annotations.NotNull;

public final class EcrfFieldReader extends EcrfReader {

    private static final String ITEM_DEF_TAG = "ItemDef";
    private static final String ITEM_DEF_OID_ATTRIBUTE = "OID";
    private static final String ITEM_DEF_NAME_ATTRIBUTE = "Name";
    private static final String ITEM_DEF_CODE_LIST_REF = "CodeListRef";
    private static final String ITEM_DEF_CODE_LIST_OID = "CodeListOID";

    private static final String CODE_LIST_TAG = "CodeList";
    private static final String CODE_LIST_OID_ATTRIBUTE = "OID";
    private static final String CODE_LIST_ITEM_TAG = "CodeListItem";
    private static final String CODE_LIST_ITEM_VALUE_ATTRIBUTE = "CodedValue";
    private static final String CODE_LIST_ITEM_CONTENT = "TranslatedText";

    private EcrfFieldReader() {
    }

    @NotNull
    public static List<EcrfField> readFields(@NotNull XMLStreamReader reader) throws XMLStreamException {
        ODMContainer container = extractODM(reader);
        return odmToEcrfFields(container);
    }

    @NotNull
    @VisibleForTesting
    static ODMContainer extractODM(@NotNull XMLStreamReader reader) throws XMLStreamException {
        List<ItemDef> itemDefs = Lists.newArrayList();
        List<CodeList> codeLists = Lists.newArrayList();
        while (reader.hasNext() && !isClinicalDataStart(reader)) {
            if (isItemDefStart(reader)) {
                itemDefs.add(extractItemDef(reader));
            }
            if (isCodeListStart(reader)) {
                codeLists.add(extractCodeList(reader));
            }
            next(reader);
        }

        return new ODMContainer(itemDefs, codeLists);
    }

    @NotNull
    @VisibleForTesting
    static List<EcrfField> odmToEcrfFields(@NotNull ODMContainer container) {
        List<EcrfField> fields = Lists.newArrayListWithCapacity(container.itemDefs().size());
        for (ItemDef item : container.itemDefs()) {
            Map<Integer, String> values = Maps.newHashMap();
            if (item.codeListOID() != null) {
                values = findValuesForCodeList(container.codeLists(), item.codeListOID());
            }
            String category = OIDFunctions.category(item.OID());
            String fieldName = OIDFunctions.fieldName(item.OID());
            String description = item.name();

            fields.add(new EcrfField(category, fieldName, description, values));
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

    @NotNull
    private static ItemDef extractItemDef(@NotNull XMLStreamReader reader) throws XMLStreamException {
        String OID = reader.getAttributeValue("", ITEM_DEF_OID_ATTRIBUTE);
        String name = reader.getAttributeValue("", ITEM_DEF_NAME_ATTRIBUTE);
        String codeListOID = null;
        while (!isItemDefEnd(reader)) {
            next(reader);
            if (isCodeListRefStart(reader)) {
                codeListOID = reader.getAttributeValue("", ITEM_DEF_CODE_LIST_OID);
            }
        }

        return new ItemDef(OID, name, codeListOID);
    }

    @NotNull
    private static CodeList extractCodeList(@NotNull XMLStreamReader reader) throws XMLStreamException {
        String OID = reader.getAttributeValue("", CODE_LIST_OID_ATTRIBUTE);
        Map<Integer, String> codeListItems = Maps.newHashMap();
        while (!isCodeListEnd(reader)) {
            next(reader);
            if (isCodeListItem(reader)) {
                int index = Integer.valueOf(reader.getAttributeValue("", CODE_LIST_ITEM_VALUE_ATTRIBUTE));
                while (!isCodeListItemContent(reader)) {
                    next(reader);
                }
                // KODU: The item content for the code lists is held in text in the next event!
                next(reader);
                codeListItems.put(index, CodeListFactory.fromText(reader.getText()));
            }
        }
        return new CodeList(OID, codeListItems);
    }

    private static boolean isItemDefStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_DEF_TAG);
    }

    private static boolean isItemDefEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, ITEM_DEF_TAG);
    }

    private static boolean isCodeListRefStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, ITEM_DEF_CODE_LIST_REF);
    }

    private static boolean isCodeListStart(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListEnd(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.END_ELEMENT, CODE_LIST_TAG);
    }

    private static boolean isCodeListItem(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_ITEM_TAG);
    }

    private static boolean isCodeListItemContent(@NotNull XMLStreamReader reader) {
        return isOfTypeWithName(reader, XMLEvent.START_ELEMENT, CODE_LIST_ITEM_CONTENT);
    }
}
