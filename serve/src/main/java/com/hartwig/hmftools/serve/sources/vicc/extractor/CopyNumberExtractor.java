package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.MutationType;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.serve.sources.vicc.check.GeneChecker;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    @NotNull
    private final Map<String, HmfTranscriptRegion> transcriptPerGeneMap;
    @NotNull
    private final GeneChecker geneChecker;

    public CopyNumberExtractor(@NotNull final Map<String, HmfTranscriptRegion> transcriptPerGeneMap,
            @NotNull final GeneChecker geneChecker) {
        this.transcriptPerGeneMap = transcriptPerGeneMap;
        this.geneChecker = geneChecker;
    }

    @NotNull
    public Map<Feature, KnownCopyNumber> extractAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownCopyNumber> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == MutationType.AMPLIFICATION || feature.type() == MutationType.DELETION) {
                HmfTranscriptRegion canonicalTranscript = transcriptPerGeneMap.get(feature.geneSymbol());

                if (geneChecker.isValidGene(feature.geneSymbol(), canonicalTranscript, feature.name(), null)) {
                    CopyNumberType type =
                            feature.type() == MutationType.AMPLIFICATION ? CopyNumberType.AMPLIFICATION : CopyNumberType.DELETION;
                    ampsDelsPerFeature.put(feature, ImmutableKnownCopyNumber.builder().gene(feature.geneSymbol()).type(type).build());
                }
            }
        }

        return ampsDelsPerFeature;
    }
}
