package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;
import java.util.Objects;
import java.util.Optional;

import com.google.common.collect.Lists;

import org.jetbrains.annotations.NotNull;

public final class XMLEcrfChecker {

    private XMLEcrfChecker() {
    }

    @NotNull
    public static List<String> checkReferences(@NotNull final XMLEcrfDatamodel datamodel) {
        final List<String> missingDefs = Lists.newArrayList();
        datamodel.studyEvents().stream().flatMap(studyEvent -> studyEvent.formOIDs().stream()).forEach(formOID -> {
            final Optional<Form> formOpt = findByOID(datamodel.forms(), formOID);
            if (!formOpt.isPresent()) {
                missingDefs.add("Missing Form: " + formOID);
            } else
                formOpt.ifPresent(form -> form.itemGroupOIDs().forEach(itemGroupOID -> {
                    final Optional<ItemGroup> itemGroupOpt = findByOID(datamodel.itemGroups(), itemGroupOID);
                    if (!itemGroupOpt.isPresent()) {
                        missingDefs.add("Missing ItemGroup: " + itemGroupOID);
                    } else
                        itemGroupOpt.ifPresent(itemGroup -> itemGroup.itemOIDs().forEach(itemOID -> {
                            final Optional<Item> itemOpt = findByOID(datamodel.items(), itemOID);
                            if (!itemOpt.isPresent()) {
                                missingDefs.add("Missing Item: " + itemOID);
                            } else
                                itemOpt.map(Item::codeListOID).filter(Objects::nonNull).ifPresent(codeListOID -> {
                                    if (!findByOID(datamodel.codeLists(), codeListOID).isPresent()) {
                                        missingDefs.add("Missing CodeList: " + codeListOID);
                                    }
                                });
                        }));
                }));
        });
        return missingDefs;
    }

    @NotNull
    private static <T extends OIDObject> Optional<T> findByOID(@NotNull final List<T> objects,
            @NotNull final String OID) {
        for (final T object : objects) {
            if (object.OID().equals(OID)) {
                return Optional.of(object);
            }
        }
        return Optional.empty();
    }
}
