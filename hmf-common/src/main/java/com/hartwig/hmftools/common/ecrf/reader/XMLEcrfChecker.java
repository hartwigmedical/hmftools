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
        datamodel.studyEvents().values().stream().flatMap(studyEvent -> studyEvent.formOIDs().stream()).forEach(formOID -> {
            final Optional<Form> formOpt = Optional.ofNullable(datamodel.forms().get(formOID));
            if (!formOpt.isPresent()) {
                missingDefs.add("Missing Form: " + formOID);
            } else {
                formOpt.ifPresent(form -> form.itemGroupOIDs().forEach(itemGroupOID -> {
                    final Optional<ItemGroup> itemGroupOpt = Optional.ofNullable(datamodel.itemGroups().get(itemGroupOID));
                    if (!itemGroupOpt.isPresent()) {
                        missingDefs.add("Missing ItemGroup: " + itemGroupOID);
                    } else {
                        itemGroupOpt.ifPresent(itemGroup -> itemGroup.itemOIDs().forEach(itemOID -> {
                            final Optional<Item> itemOpt = Optional.ofNullable(datamodel.items().get(itemOID));
                            if (!itemOpt.isPresent()) {
                                missingDefs.add("Missing Item: " + itemOID);
                            } else {
                                itemOpt.map(Item::codeListOID).filter(Objects::nonNull).ifPresent(codeListOID -> {
                                    if (!Optional.ofNullable(datamodel.codeLists().get(codeListOID)).isPresent()) {
                                        missingDefs.add("Missing CodeList: " + codeListOID);
                                    }
                                });
                            }
                        }));
                    }
                }));
            }
        });
        return missingDefs;
    }
}
