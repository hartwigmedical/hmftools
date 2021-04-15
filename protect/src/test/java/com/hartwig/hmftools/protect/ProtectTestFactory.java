package com.hartwig.hmftools.protect;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.protect.ImmutableProtectEvidence;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.actionability.EvidenceDirection;
import com.hartwig.hmftools.protect.evidence.PersonalizedEvidenceFactory;
import com.hartwig.hmftools.protect.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.germline.ImmutableGermlineReportingEntry;
import com.hartwig.hmftools.serve.actionability.ActionabilityTestUtil;
import com.hartwig.hmftools.serve.actionability.ActionableEvent;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProtectTestFactory {

    private ProtectTestFactory() {
    }

    @NotNull
    public static PersonalizedEvidenceFactory createTestEvidenceFactory() {
        return new PersonalizedEvidenceFactory(Sets.newHashSet());
    }

    @NotNull
    public static ActionableEvent createTestEvent() {
        return ActionabilityTestUtil.create(Knowledgebase.VICC_CGI,
                Strings.EMPTY,
                Strings.EMPTY,
                Strings.EMPTY,
                com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A,
                EvidenceDirection.RESPONSIVE,
                Sets.newHashSet());
    }

    @NotNull
    public static ImmutableProtectEvidence.Builder createTestBuilder() {
        return ImmutableProtectEvidence.builder()
                .genomicEvent(Strings.EMPTY)
                .germline(false)
                .reported(true)
                .treatment(Strings.EMPTY)
                .onLabel(false)
                .level(com.hartwig.hmftools.common.serve.actionability.EvidenceLevel.A)
                .direction(EvidenceDirection.RESPONSIVE);
    }

    @NotNull
    public static GermlineReportingModel createTestGermlineModel(@NotNull String gene, boolean reportBiallelicOnly,
            @Nullable String exclusiveHgvsProteinFilter) {
        return new GermlineReportingModel(Lists.newArrayList(ImmutableGermlineReportingEntry.builder()
                .gene(gene)
                .notifyClinicalGeneticist(false)
                .exclusiveHgvsProteinFilter(exclusiveHgvsProteinFilter)
                .build()));
    }

    @NotNull
    public static GermlineReportingModel createEmptyGermlineReportingModel() {
        return new GermlineReportingModel(Lists.newArrayList());
    }
}
