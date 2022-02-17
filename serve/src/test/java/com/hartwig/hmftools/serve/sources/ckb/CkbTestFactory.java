package com.hartwig.hmftools.serve.sources.ckb;

import java.time.LocalDate;

import com.hartwig.hmftools.ckb.datamodel.indication.ImmutableIndication;
import com.hartwig.hmftools.ckb.datamodel.indication.Indication;
import com.hartwig.hmftools.ckb.datamodel.therapy.Therapy;
import com.hartwig.hmftools.ckb.datamodel.CkbEntry;
import com.hartwig.hmftools.ckb.datamodel.ImmutableCkbEntry;
import com.hartwig.hmftools.ckb.datamodel.evidence.Evidence;
import com.hartwig.hmftools.ckb.datamodel.evidence.ImmutableEvidence;
import com.hartwig.hmftools.ckb.datamodel.therapy.ImmutableTherapy;
import com.hartwig.hmftools.ckb.datamodel.variant.Gene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableGene;
import com.hartwig.hmftools.ckb.datamodel.variant.ImmutableVariant;
import com.hartwig.hmftools.ckb.datamodel.variant.Variant;

import org.apache.commons.compress.utils.Lists;
import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;

public final class CkbTestFactory {

    private static final LocalDate TEST_DATE = LocalDate.of(2021, 2, 20);

    private CkbTestFactory() {
    }

    @NotNull
    public static CkbEntry createEntryWithGene(@NotNull String geneSymbol) {
        return createEntry(geneSymbol,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntryWithVariant(@NotNull String variant) {
        return createEntry(Strings.EMPTY,
                variant,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntryWithFullName(@NotNull String fullName) {
        return createEntry(Strings.EMPTY,
                Strings.EMPTY,
                fullName,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntryWithGeneAndVariant(@NotNull String geneSymbol, @NotNull String variant) {
        return createEntry(geneSymbol,
                variant,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY);
    }

    @NotNull
    public static CkbEntry createEntry(@NotNull String geneSymbol, @NotNull String variant, @NotNull String fullName,
            @NotNull String responseType, @NotNull String evidenceType, @NotNull String therapyName, @NotNull String name,
            @NotNull String level, @NotNull String termId) {
        return ImmutableCkbEntry.builder()
                .profileId(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .profileName(Strings.EMPTY)
                .addVariants(createVariant(geneSymbol, variant, fullName))
                .addEvidences(createEvidence(responseType, evidenceType, therapyName, name, level, termId))
                .build();
    }

    @NotNull
    private static Evidence createEvidence(@NotNull String responseType, @NotNull String evidenceType, @NotNull String therapyName,
            @NotNull String name, @NotNull String level, @NotNull String termId) {
        return ImmutableEvidence.builder()
                .id(0)
                .therapy(createTherapy(therapyName))
                .indication(createIndication(name, termId))
                .responseType(responseType)
                .relevantTreatmentApproaches(Lists.newArrayList())
                .evidenceType(evidenceType)
                .efficacyEvidence(Strings.EMPTY)
                .approvalStatus(Strings.EMPTY)
                .ampCapAscoEvidenceLevel(level)
                .ampCapAscoInferredTier(Strings.EMPTY)
                .references(Lists.newArrayList())
                .build();
    }

    @NotNull
    private static Therapy createTherapy(@NotNull String therapyName) {
        return ImmutableTherapy.builder()
                .id(0)
                .createDate(TEST_DATE)
                .updateDate(TEST_DATE)
                .therapyName(therapyName)
                .drugs(Lists.newArrayList())
                .synonyms(Lists.newArrayList())
                .globalApprovalStatuses(Lists.newArrayList())
                .description(Strings.EMPTY)
                .references(Lists.newArrayList())
                .build();
    }

    @NotNull
    private static Indication createIndication(@NotNull String name, @NotNull String termId) {
        return ImmutableIndication.builder()
                .id(0)
                .name(name)
                .source(Strings.EMPTY)
                .definition(Strings.EMPTY)
                .currentPreferredTerm(Strings.EMPTY)
                .lastUpdateDateFromDO(TEST_DATE)
                .termId(termId)
                .altIds(Lists.newArrayList())
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