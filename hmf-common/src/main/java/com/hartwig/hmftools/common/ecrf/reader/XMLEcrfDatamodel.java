package com.hartwig.hmftools.common.ecrf.reader;

import java.util.List;

import org.immutables.value.Value;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

@Value.Immutable
@Value.Style(of = "new",
             allParameters = true,
             passAnnotations = { NotNull.class, Nullable.class })
public abstract class XMLEcrfDatamodel {

    @NotNull
    public abstract List<StudyEvent> studyEvents();

    @NotNull
    public abstract List<Form> forms();

    @NotNull
    public abstract List<ItemGroup> itemGroups();

    @NotNull
    public abstract List<Item> items();

    @NotNull
    public abstract List<CodeList> codeLists();
}
