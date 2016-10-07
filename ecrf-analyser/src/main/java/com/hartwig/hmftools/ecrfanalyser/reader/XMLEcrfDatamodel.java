package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

public class XMLEcrfDatamodel {

    @NotNull
    private final List<StudyEvent> studyEvents;
    @NotNull
    private final List<Form> forms;
    @NotNull
    private final List<ItemGroup> itemGroups;
    @NotNull
    private final List<Item> items;
    @NotNull
    private final List<CodeList> codeLists;

    XMLEcrfDatamodel(@NotNull final List<StudyEvent> studyEvents, @NotNull final List<Form> forms,
            @NotNull final List<ItemGroup> itemGroups, @NotNull final List<Item> items,
            @NotNull final List<CodeList> codeLists) {
        this.studyEvents = studyEvents;
        this.forms = forms;
        this.itemGroups = itemGroups;
        this.items = items;
        this.codeLists = codeLists;
    }

    @NotNull
    public List<StudyEvent> studyEvents() {
        return studyEvents;
    }

    @NotNull
    public List<Form> forms() {
        return forms;
    }

    @NotNull
    public List<ItemGroup> itemGroups() {
        return itemGroups;
    }

    @NotNull
    List<Item> items() {
        return items;
    }

    @NotNull
    List<CodeList> codeLists() {
        return codeLists;
    }
}
