package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.ViccUtil;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    private static final Set<EventType> COPY_NUMBER_MUTATIONS = Sets.newHashSet(EventType.AMPLIFICATION, EventType.DELETION);

    @NotNull
    private final GeneChecker geneChecker;

    public CopyNumberExtractor(@NotNull final GeneChecker geneChecker) {
        this.geneChecker = geneChecker;
    }

    @NotNull
    public Map<Feature, KnownCopyNumber> extract(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownCopyNumber> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (COPY_NUMBER_MUTATIONS.contains(feature.type()) && geneChecker.isValidGene(feature.geneSymbol())) {
                CopyNumberType type = feature.type() == EventType.AMPLIFICATION ? CopyNumberType.AMPLIFICATION : CopyNumberType.DELETION;

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
