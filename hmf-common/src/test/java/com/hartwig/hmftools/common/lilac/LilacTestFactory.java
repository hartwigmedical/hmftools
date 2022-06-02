package com.hartwig.hmftools.common.lilac;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LilacTestFactory {

    private LilacTestFactory() {
    }

    @NotNull
    public static ImmutableLilacRecord.Builder builder() {
        return ImmutableLilacRecord.builder()
                .allele(Strings.EMPTY)
                .refFragments(0)
                .tumorFragments(0)
                .tumorCopyNumber(0D)
                .somaticMissense(0D)
                .somaticNonsenseOrFrameshift(0D)
                .somaticSplice(0D)
                .somaticSynonymous(0D)
                .somaticInframeIndel(0D);
    }
}
