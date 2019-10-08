package com.hartwig.hmftools.common.ecrf.reader;

import static java.util.stream.Collectors.toMap;

import java.util.List;
import java.util.Map;
import java.util.function.Function;

import com.hartwig.hmftools.common.ecrf.datamodel.EcrfResolveException;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class XMLEcrfDatamodel {
    private static final String BIRTH_DATE_IDENTIFIER = "BIRTHDTC";

    @NotNull
    public abstract Map<String, StudyEvent> studyEvents();

    @NotNull
    public abstract Map<String, Form> forms();

    @NotNull
    public abstract Map<String, ItemGroup> itemGroups();

    @NotNull
    public abstract Map<String, Item> items();

    @NotNull
    public abstract Map<String, CodeList> codeLists();

    @NotNull
    public static XMLEcrfDatamodel of(@NotNull final List<StudyEvent> studyEvents, @NotNull final List<Form> forms,
            @NotNull final List<ItemGroup> itemGroups, @NotNull final List<Item> items, @NotNull final List<CodeList> codeLists) {
        final Map<String, StudyEvent> studyEventMap = studyEvents.stream().collect(toMap(OIDObject::oid, Function.identity()));
        final Map<String, Form> formMap = forms.stream().collect(toMap(OIDObject::oid, Function.identity()));
        final Map<String, ItemGroup> itemGroupMap = itemGroups.stream().collect(toMap(OIDObject::oid, Function.identity()));
        final Map<String, Item> itemMap = items.stream().collect(toMap(OIDObject::oid, Function.identity()));
        final Map<String, CodeList> codeListMap = codeLists.stream().collect(toMap(OIDObject::oid, Function.identity()));
        return ImmutableXMLEcrfDatamodel.of(studyEventMap, formMap, itemGroupMap, itemMap, codeListMap);
    }

    @NotNull
    String resolveValue(@NotNull final String itemOID, @NotNull final String ecrfValue) throws EcrfResolveException {
        String value;
        if (items().get(itemOID) != null && items().get(itemOID).codeListOID() != null) {
            final String codeListOID = items().get(itemOID).codeListOID();
            final CodeList codeList = codeLists().get(codeListOID);
            if (codeList != null && codeList.values().size() > 0 && ecrfValue.length() > 0) {
                if (isInteger(ecrfValue)) {
                    value = codeList.values().get(Integer.valueOf(ecrfValue));
                    if (value == null) {
                        throw new EcrfResolveException("Could not find value in dropdown for " + itemOID + ": " + ecrfValue);
                    }
                } else {
                    throw new EcrfResolveException("Could not convert value from list to integer for " + itemOID + ": " + ecrfValue);
                }
            } else {
                value = ecrfValue;
            }

            if (itemOID.contains(BIRTH_DATE_IDENTIFIER) && value.length() > 0) {
                if (value.length() < 4) {
                    throw new EcrfResolveException("Could not convert " + itemOID + ": " + value);
                } else {
                    value = value.substring(0, 4) + "-01-01";
                }
            }
            return value;
        }

        return ecrfValue;
    }

    private static boolean isInteger(@NotNull String integerString) {
        try {
            Integer.valueOf(integerString);
            return true;
        } catch (NumberFormatException exception) {
            return false;
        }
    }
}
