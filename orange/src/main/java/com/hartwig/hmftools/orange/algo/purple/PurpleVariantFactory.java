package com.hartwig.hmftools.orange.algo.purple;

import java.util.List;

import com.hartwig.hmftools.common.variant.SomaticVariant;

import org.apache.commons.compress.utils.Lists;
import org.jetbrains.annotations.Nullable;

public final class PurpleVariantFactory {

    private PurpleVariantFactory() {
    }

    @Nullable
    public static List<PurpleVariant> create(@Nullable List<SomaticVariant> variants) {
        if (variants == null) {
            return null;
        }

        // TODO Implement
        return Lists.newArrayList();
    }

//    @NotNull
//    public static List<PurpleVariant> merge(@NotNull List<PurpleVariant> variants1, @NotNull List<PurpleVariant> variants2) {
//        List<PurpleVariant> merged = Lists.newArrayList();
//        merged.addAll(variants1);
//        merged.addAll(variants2);
//        return merged;
//    }
}
