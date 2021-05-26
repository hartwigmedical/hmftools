package com.hartwig.hmftools.serve.sources.ckb;

import java.time.LocalDate;

import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.variant.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class CkbTestFactory {

    private static final LocalDate TEST_DATE = LocalDate.of(2021, 2, 20);

    private CkbTestFactory() {
    }

    @NotNull
    public static CkbEntry createEntryWithGene(@NotNull String geneSymbol) {
        return createEntry(geneSymbol, Strings.EMPTY, Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntryWithVariant(@NotNull String variant) {
        return createEntry(Strings.EMPTY, variant, Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntryWithFullName(@NotNull String fullName) {
        return createEntry(Strings.EMPTY, Strings.EMPTY, fullName);
    }

    @NotNull
    public static CkbEntry createEntryWithGeneAndVariant(@NotNull String geneSymbol, @NotNull String variant) {
        return createEntry(geneSymbol, variant, Strings.EMPTY);
    }

    @NotNull
    private static CkbEntry createEntry(@NotNull String geneSymbol, @NotNull String variant, @NotNull String fullName) {
        return ImmutableCkbEntry.builder()
                .profileId(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .profileName(Strings.EMPTY)
                .addVariants(createVariant(geneSymbol, variant, fullName))
                .build();
    }

    @NotNull
    private static Variant createVariant(@NotNull String geneSymbol, @NotNull String variant, @NotNull String fullName) {
        return ImmutableVariant.builder()
                .id(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .fullName(fullName)
                .variant(variant)
                .gene(createGene(geneSymbol))
                .build();
    }

    @NotNull
    private static Gene createGene(@NotNull String geneSymbol) {
        return ImmutableGene.builder()
                .id(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .geneSymbol(geneSymbol)
                .geneRole(Strings.EMPTY)
                .build();
    }
}
