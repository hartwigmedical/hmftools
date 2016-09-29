package com.hartwig.hmftools.ecrfanalyser.reader;

import org.jetbrains.annotations.NotNull;

final class ItemDefTestFunctions {

    @NotNull
    static String toOID(@NotNull String category, @NotNull String name) {
        return "FLD." + category + "." + name;
    }
}
