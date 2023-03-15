package com.hartwig.hmftools.datamodel.lilac;

import com.hartwig.hmftools.datamodel.hla.ImmutableLilacAllele;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class TestLilacFactory {

    private TestLilacFactory() {
    }

    @NotNull
    public static ImmutableLilacAllele.Builder builder() {
        return ImmutableLilacAllele.builder()
                .allele(Strings.EMPTY)
                .tumorCopyNumber(0D)
                .refFragments(0)
                .tumorFragments(0)
                .rnaFragments(0)
                .somaticMissense(0D)
                .somaticNonsenseOrFrameshift(0D)
                .somaticSplice(0D)
                .somaticSynonymous(0D)
                .somaticInframeIndel(0D);
    }
}
