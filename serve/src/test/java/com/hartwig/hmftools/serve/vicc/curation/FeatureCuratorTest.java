package com.hartwig.hmftools.serve.vicc.curation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.serve.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class FeatureCuratorTest {

    @Test
    public void canCurateFeatures() {
        FeatureCurator.CurationKey firstKey = FeatureCurator.ONCOKB_FEATURE_NAME_MAPPINGS.keySet().iterator().next();
        FeatureNameMapping firstMapping = FeatureCurator.ONCOKB_FEATURE_NAME_MAPPINGS.get(firstKey).get(0);

        ViccEntry entry =
                ViccTestFactory.testViccEntryForSource(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript(firstKey.transcript()));

        Feature feature = ImmutableFeature.builder().geneSymbol(firstKey.gene()).name(firstMapping.originalFeatureName()).build();

        assertEquals(firstMapping.curatedFeatureName(), FeatureCurator.curate(entry, feature).name());
    }
}