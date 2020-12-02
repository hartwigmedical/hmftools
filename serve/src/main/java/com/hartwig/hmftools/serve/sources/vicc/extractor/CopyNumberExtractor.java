package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.ViccUtil;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    private static final Set<MutationType> COPY_NUMBER_MUTATIONS = Sets.newHashSet(MutationType.AMPLIFICATION, MutationType.DELETION);

    @NotNull
    private final GeneChecker geneChecker;

    public CopyNumberExtractor(@NotNull final GeneChecker geneChecker) {
        this.geneChecker = geneChecker;
    }

    @NotNull
    public Map<Feature, KnownCopyNumber> extractAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownCopyNumber> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (geneChecker.isValidGene(feature.geneSymbol()) && COPY_NUMBER_MUTATIONS.contains(feature.type())) {
                CopyNumberType type = feature.type() == MutationType.AMPLIFICATION ? CopyNumberType.AMPLIFICATION : CopyNumberType.DELETION;

                ampsDelsPerFeature.put(feature,
                        ImmutableKnownCopyNumber.builder()
                                .gene(feature.geneSymbol())
                                .type(type)
                                .addSources(ViccUtil.toKnowledgebase(viccEntry.source()))
                                .build());

            }
        }

        return ampsDelsPerFeature;
    }
}
