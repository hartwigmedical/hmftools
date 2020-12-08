package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Set;

import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.Knowledgebase;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public class CopyNumberExtractor {

    private static final Set<EventType> COPY_NUMBER_MUTATIONS = Sets.newHashSet(EventType.AMPLIFICATION, EventType.DELETION);

    @NotNull
    private final GeneChecker geneChecker;

    public CopyNumberExtractor(@NotNull final GeneChecker geneChecker) {
        this.geneChecker = geneChecker;
    }

    @Nullable
    public KnownCopyNumber extract(@NotNull Knowledgebase source, @NotNull String gene, @NotNull EventType type) {
        if (COPY_NUMBER_MUTATIONS.contains(type) && geneChecker.isValidGene(gene)) {
            CopyNumberType copyNumberType = type == EventType.AMPLIFICATION ? CopyNumberType.AMPLIFICATION : CopyNumberType.DELETION;

            return ImmutableKnownCopyNumber.builder().gene(gene).type(copyNumberType).addSources(source).build();
        }

        return null;
    }
}
