package com.hartwig.hmftools.serve.sources.vicc.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class FeatureCuratorTest {

    @Test
    public void canCurateFeatures() {
        FeatureCurationKey firstOncoKbKey = firstOncoKbMappingKey();
        String firstMappedFeature = FeatureCurationFactory.FEATURE_MAPPINGS.get(firstOncoKbKey).featureName();

        ViccEntry entry = ViccTestFactory.testEntryWithSourceAndTranscript(ViccSource.ONCOKB, firstOncoKbKey.transcript());

        Feature feature = ViccTestFactory.testFeatureWithGeneAndName(firstOncoKbKey.gene(), firstOncoKbKey.featureName());

        assertEquals(firstMappedFeature, new FeatureCurator().curate(entry, feature).name());
    }

    @Test
    public void canBlacklistFeatures() {
        FeatureCurationKey firstOncoKbKey = firstOncoKbBlacklistKey();
        ViccEntry entry = ViccTestFactory.testEntryWithSourceAndTranscript(ViccSource.ONCOKB, firstOncoKbKey.transcript());

        Feature feature = ViccTestFactory.testFeatureWithGeneAndName(firstOncoKbKey.gene(), firstOncoKbKey.featureName());
        assertNull(new FeatureCurator().curate(entry, feature));
    }

    @Test
    public void canKeepTrackOfFeatures() {
        FeatureCurator curator = new FeatureCurator();

        ViccEntry entry = ViccTestFactory.testEntryWithSourceAndTranscript(ViccSource.ONCOKB, "any");
        Feature feature = ViccTestFactory.testFeatureWithGeneAndName("any", "any");

        assertNotNull(curator.curate(entry, feature));

        FeatureCurationKey blacklistKey = firstOncoKbBlacklistKey();
        ViccEntry blacklistEntry = ViccTestFactory.testEntryWithSourceAndTranscript(ViccSource.ONCOKB, blacklistKey.transcript());

        Feature blacklistFeature = ViccTestFactory.testFeatureWithGeneAndName(blacklistKey.gene(), blacklistKey.featureName());

        assertNull(curator.curate(blacklistEntry, blacklistFeature));

        curator.reportUnusedCurationKeys();
    }

    @NotNull
    private static FeatureCurationKey firstOncoKbMappingKey() {
        for (FeatureCurationKey key : FeatureCurationFactory.FEATURE_MAPPINGS.keySet()) {
            if (key.source() == ViccSource.ONCOKB) {
                return key;
            }
        }
        throw new IllegalStateException("No OncoKB mapping keys found!");
    }

    @NotNull
    private static FeatureCurationKey firstOncoKbBlacklistKey() {
        for (FeatureCurationKey key : FeatureCurationFactory.FEATURE_BLACKLIST) {
            if (key.source() == ViccSource.ONCOKB) {
                return key;
            }
        }
        throw new IllegalStateException("No OncoKB blacklist keys found!");
    }
}