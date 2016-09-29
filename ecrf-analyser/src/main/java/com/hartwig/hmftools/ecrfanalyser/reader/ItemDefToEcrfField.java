package com.hartwig.hmftools.ecrfanalyser.reader;

import org.jetbrains.annotations.NotNull;

final class ItemDefToEcrfField {

    private static final String ITEM_DEF_OID_SEPARATOR = "\\.";
    private static final String NO_CATEGORY = "GENERAL";

    private ItemDefToEcrfField() {
    }

    @NotNull
    static String category(@NotNull ItemDef itemDef) {
        String[] fields = itemDef.OID().split(ITEM_DEF_OID_SEPARATOR);
        if (fields.length == 2) {
            return NO_CATEGORY;
        }

        assert fields.length == 3;
        return fields[1];
    }

    @NotNull
    static String fieldName(@NotNull ItemDef itemDef) {
        String[] fields = itemDef.OID().split(ITEM_DEF_OID_SEPARATOR);
        return fields[fields.length - 1];
    }

    @NotNull
    static String description(@NotNull ItemDef itemDef) {
        return itemDef.name();
    }
}
