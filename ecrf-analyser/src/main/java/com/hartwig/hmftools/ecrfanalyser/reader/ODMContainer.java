package com.hartwig.hmftools.ecrfanalyser.reader;

import java.util.List;

import org.jetbrains.annotations.NotNull;

class ODMContainer {

    @NotNull
    private final List<ItemDef> itemDefs;
    @NotNull
    private final List<CodeList> codeLists;

    ODMContainer(@NotNull final List<ItemDef> itemDefs, @NotNull final List<CodeList> codeLists) {
        this.itemDefs = itemDefs;
        this.codeLists = codeLists;
    }

    @NotNull
    List<ItemDef> itemDefs() {
        return itemDefs;
    }

    @NotNull
    List<CodeList> codeLists() {
        return codeLists;
    }
}
