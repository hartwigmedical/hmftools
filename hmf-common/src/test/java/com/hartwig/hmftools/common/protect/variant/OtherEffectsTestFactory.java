package com.hartwig.hmftools.common.protect.variant;

import org.jetbrains.annotations.NotNull;

public final class OtherEffectsTestFactory {

    private OtherEffectsTestFactory() {
    }

    @NotNull
    public static String create() {
        return "trans|coding|impact|effect|MISSENSE";
    }
}
