package com.hartwig.hmftools.iclusion.classification;

import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.serve.classification.ImmutableEventClassifierConfig;

import org.jetbrains.annotations.NotNull;

public final class ClassificationConfig {

    private static final String EXON_KEYWORD = "exon";
    private static final Set<String> EXON_RANGE_EVENTS = Sets.newHashSet();
    private static final Set<String> EXON_RANGE_KEYWORDS = Sets.newHashSet();
    private static final Map<String, Set<String>> FUSION_PAIR_AND_EXON_RANGES_PER_GENE = Maps.newHashMap();
    private static final Set<String> GENERIC_GENE_LEVEL_KEYWORDS = Sets.newHashSet();
    private static final Set<String> ACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet();
    private static final Set<String> INACTIVATING_GENE_LEVEL_KEYWORDS = Sets.newHashSet();
    private static final Set<String> AMPLIFICATION_KEYWORDS = Sets.newHashSet();
    private static final Set<String> AMPLIFICATION_KEY_PHRASES = Sets.newHashSet();
    private static final Set<String> DELETION_KEYWORDS = Sets.newHashSet();
    private static final Set<String> DELETION_KEY_PHRASES = Sets.newHashSet();
    private static final Set<String> DELETION_KEYWORDS_TO_SKIP = Sets.newHashSet();
    private static final Set<String> EXONIC_DEL_DUP_FUSION_EVENTS = Sets.newHashSet();
    private static final Set<String> FUSION_PAIR_EVENTS_TO_SKIP = Sets.newHashSet();
    private static final Set<String> PROMISCUOUS_FUSION_KEYWORDS = Sets.newHashSet();
    private static final Set<String> SIGNATURE_EVENTS = Sets.newHashSet();
    private static final Map<String, Set<String>> COMBINED_EVENTS_PER_GENE = Maps.newHashMap();
    private static final Map<String, Set<String>> COMPLEX_EVENTS_PER_GENE = Maps.newHashMap();

    private ClassificationConfig() {
    }

    @NotNull
    static EventClassifierConfig build() {
        return ImmutableEventClassifierConfig.builder()
                .proteinAnnotationExtractor(new ProteinAnnotationExtractor())
                .exonKeyword(EXON_KEYWORD)
                .exonRangeEvents(EXON_RANGE_EVENTS)
                .exonRangeKeywords(EXON_RANGE_KEYWORDS)
                .fusionPairAndExonRangesPerGene(FUSION_PAIR_AND_EXON_RANGES_PER_GENE)
                .genericGeneLevelKeywords(GENERIC_GENE_LEVEL_KEYWORDS)
                .activatingGeneLevelKeywords(ACTIVATING_GENE_LEVEL_KEYWORDS)
                .inactivatingGeneLevelKeywords(INACTIVATING_GENE_LEVEL_KEYWORDS)
                .amplificationKeywords(AMPLIFICATION_KEYWORDS)
                .amplificationKeyPhrases(AMPLIFICATION_KEY_PHRASES)
                .deletionKeywords(DELETION_KEYWORDS)
                .deletionKeyPhrases(DELETION_KEY_PHRASES)
                .deletionKeywordsToSkip(DELETION_KEYWORDS_TO_SKIP)
                .exonicDelDupFusionEvents(EXONIC_DEL_DUP_FUSION_EVENTS)
                .fusionPairEventsToSkip(FUSION_PAIR_EVENTS_TO_SKIP)
                .promiscuousFusionKeywords(PROMISCUOUS_FUSION_KEYWORDS)
                .signatureEvents(SIGNATURE_EVENTS)
                .combinedEventsPerGene(COMBINED_EVENTS_PER_GENE)
                .complexEventsPerGene(COMPLEX_EVENTS_PER_GENE)
                .build();
    }
}
