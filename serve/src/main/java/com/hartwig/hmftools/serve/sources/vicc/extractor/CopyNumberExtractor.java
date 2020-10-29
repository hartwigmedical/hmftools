package com.hartwig.hmftools.serve.sources.vicc.extractor;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.serve.copynumber.CopyNumberType;
import com.hartwig.hmftools.serve.copynumber.ImmutableKnownCopyNumber;
import com.hartwig.hmftools.serve.copynumber.KnownCopyNumber;
import com.hartwig.hmftools.vicc.annotation.FeatureType;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;

public class CopyNumberExtractor {

    public CopyNumberExtractor() {
    }

    @NotNull
    public Map<Feature, KnownCopyNumber> extractKnownAmplificationsDeletions(@NotNull ViccEntry viccEntry) {
        Map<Feature, KnownCopyNumber> ampsDelsPerFeature = Maps.newHashMap();

        for (Feature feature : viccEntry.features()) {
            if (feature.type() == FeatureType.AMPLIFICATION) {
                ampsDelsPerFeature.put(feature,
                        ImmutableKnownCopyNumber.builder().gene(feature.geneSymbol()).type(CopyNumberType.AMPLIFICATION).build());
            } else if (feature.type() == FeatureType.DELETION) {
                ampsDelsPerFeature.put(feature,
                        ImmutableKnownCopyNumber.builder().gene(feature.geneSymbol()).type(CopyNumberType.DELETION).build());
            }
        }

        return ampsDelsPerFeature;
    }
}
