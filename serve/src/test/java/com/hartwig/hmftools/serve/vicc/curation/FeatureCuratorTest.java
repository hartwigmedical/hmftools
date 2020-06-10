package com.hartwig.hmftools.serve.vicc.curation;

import static org.junit.Assert.assertEquals;

import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class FeatureCuratorTest {

    @Test
    public void canCurateFeatures() {
        String firstGene = FeatureCurator.ONCOKB_FEATURE_NAME_MAPPINGS_PER_GENE.keySet().iterator().next();
        FeatureNameMapping firstMapping = FeatureCurator.ONCOKB_FEATURE_NAME_MAPPINGS_PER_GENE.get(firstGene).get(0);

        Feature feature = ImmutableFeature.builder().geneSymbol(firstGene).name(firstMapping.originalFeatureName()).build();

        assertEquals(firstMapping.curatedFeatureName(), FeatureCurator.curate(ViccSource.ONCOKB, feature).name());
    }
}