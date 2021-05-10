package com.hartwig.hmftools.ckb;

import java.time.LocalDate;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class CkbTestFactory {

    private static final LocalDate TEST_DATE = LocalDate.of(2021, 2, 20);

    private CkbTestFactory() {
    }

    @NotNull
    public static CkbEntry createEntry(@NotNull Variant variant) {
        return createEntry(Lists.newArrayList(variant));
    }

    @NotNull
    public static CkbEntry createEntry(@NotNull List<Variant> variants) {
        return ImmutableCkbEntry.builder()
                .profileId(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .profileName(Strings.EMPTY)
                .variants(variants)
                .build();
    }

    @NotNull
    public static Variant createVariant() {
        return createVariant(Strings.EMPTY, Strings.EMPTY, Strings.EMPTY, null);
    }

    @NotNull
    public static Variant createVariant(@NotNull String fullName, @NotNull String variant, @Nullable String impact) {
        return createVariant(Strings.EMPTY, fullName, variant, impact);
    }

    @NotNull
    public static Variant createVariant(@NotNull String geneSymbol, @NotNull String fullName, @NotNull String variant,
            @Nullable String impact) {
        return createVariant(geneSymbol, Lists.newArrayList(), fullName, variant, impact);
    }

    @NotNull
    public static Variant createVariant(@NotNull String geneSymbol, @NotNull List<String> geneSynonyms, @NotNull String fullName,
            @NotNull String variant, @Nullable String impact) {
        return ImmutableVariant.builder()
                .id(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .fullName(fullName)
                .variant(variant)
                .impact(impact)
                .gene(createGene(geneSymbol, geneSynonyms))
                .build();
    }

    @NotNull
    private static Gene createGene(@NotNull String geneSymbol, @NotNull List<String> synonyms) {
        return ImmutableGene.builder()
                .id(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .geneSymbol(geneSymbol)
                .geneRole(Strings.EMPTY)
                .synonyms(synonyms)
                .build();
    }
}
