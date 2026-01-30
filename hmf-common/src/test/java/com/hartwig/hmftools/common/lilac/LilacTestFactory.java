package com.hartwig.hmftools.common.lilac;

import java.util.Collections;

import com.hartwig.hmftools.common.hla.ImmutableLilacAllele;
import com.hartwig.hmftools.common.hla.ImmutableLilacSummaryData;
import com.hartwig.hmftools.common.hla.LilacSummaryData;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class LilacTestFactory
{
    @NotNull
    public static LilacSummaryData createEmptyData()
    {
        return ImmutableLilacSummaryData.builder().qc(Collections.emptyList()).build();
    }

    @NotNull
    public static ImmutableLilacAllele.Builder alleleBuilder()
    {
        return ImmutableLilacAllele.builder()
                .allele(Strings.EMPTY)
                .refFragments(0)
                .refUnique(0)
                .refShared(0)
                .refWild(0)
                .tumorFragments(0)
                .tumorUnique(0)
                .tumorShared(0)
                .tumorWild(0)
                .rnaFragments(0)
                .rnaUnique(0)
                .rnaShared(0)
                .rnaWild(0)
                .tumorCopyNumber(0D)
                .somaticMissense(0D)
                .somaticNonsenseOrFrameshift(0D)
                .somaticSplice(0D)
                .somaticSynonymous(0D)
                .somaticInframeIndel(0D);
    }
}
