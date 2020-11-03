package com.hartwig.hmftools.protect;

import java.io.IOException;

import com.google.common.collect.Lists;
import com.google.common.io.Resources;
import com.hartwig.hmftools.common.actionability.ActionabilityAnalyzer;
import com.hartwig.hmftools.common.actionability.ActionabilitySource;
import com.hartwig.hmftools.common.actionability.EvidenceLevel;
import com.hartwig.hmftools.common.actionability.EvidenceScope;
import com.hartwig.hmftools.common.actionability.ImmutableEvidenceItem;
import com.hartwig.hmftools.protect.variants.germline.GermlineReportingModel;
import com.hartwig.hmftools.protect.variants.germline.ImmutableGermlineReportingEntry;

import org.apache.logging.log4j.util.Strings;
import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class ProtectTestFactory {

    private static final String KNOWLEDGEBASE_DIRECTORY = Resources.getResource("actionability").getPath();

    private ProtectTestFactory() {
    }

    @NotNull
    public static ActionabilityAnalyzer loadTestActionabilityAnalyzer() {
        try {
            return ActionabilityAnalyzer.fromKnowledgebase(KNOWLEDGEBASE_DIRECTORY);
        } catch (IOException e) {
            throw new IllegalStateException(e);
        }
    }

    @NotNull
    public static ImmutableEvidenceItem.Builder createTestEvidenceBuilder() {
        return ImmutableEvidenceItem.builder()
                .event(Strings.EMPTY)
                .source(ActionabilitySource.CIVIC)
                .reference(Strings.EMPTY)
                .drug(Strings.EMPTY)
                .drugsType(Strings.EMPTY)
                .level(EvidenceLevel.LEVEL_A)
                .response(Strings.EMPTY)
                .isOnLabel(false)
                .cancerType(Strings.EMPTY)
                .scope(EvidenceScope.SPECIFIC);
    }

    @NotNull
    public static GermlineReportingModel createTestGermlineModel(@NotNull String gene, boolean reportBiallelicOnly,
            @Nullable String exclusiveHgvsProteinFilter) {
        return new GermlineReportingModel(Lists.newArrayList(ImmutableGermlineReportingEntry.builder()
                .gene(gene)
                .notifyClinicalGeneticist(false)
                .reportBiallelicOnly(reportBiallelicOnly)
                .exclusiveHgvsProteinFilter(exclusiveHgvsProteinFilter)
                .build()));
    }

    @NotNull
    public static GermlineReportingModel createEmptyGermlineReportingModel() {
        return new GermlineReportingModel(Lists.newArrayList());
    }
}
