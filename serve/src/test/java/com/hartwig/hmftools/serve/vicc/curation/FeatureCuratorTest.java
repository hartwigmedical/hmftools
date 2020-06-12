package com.hartwig.hmftools.serve.vicc.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.serve.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.junit.Test;

public class FeatureCuratorTest {

    @Test
    public void canCurateFeatures() {
        CurationKey firstKey = CurationFactory.ONCOKB_FEATURE_NAME_MAPPINGS.keySet().iterator().next();
        String firstMappedFeature = CurationFactory.ONCOKB_FEATURE_NAME_MAPPINGS.get(firstKey);

        ViccEntry entry =
                ViccTestFactory.testViccEntryForSource(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript(firstKey.transcript()));

        Feature feature = ImmutableFeature.builder().geneSymbol(firstKey.gene()).name(firstKey.featureName()).build();

        assertEquals(firstMappedFeature, new FeatureCurator().curate(entry, feature).name());
    }

    @Test
    public void canBlacklistFeatures() {
        CurationKey firstKey = CurationFactory.ONCOKB_FEATURE_BLACKLIST.iterator().next();
        ViccEntry entry =
                ViccTestFactory.testViccEntryForSource(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript(firstKey.transcript()));

        Feature feature = ImmutableFeature.builder().geneSymbol(firstKey.gene()).name(firstKey.featureName()).build();
        assertNull(new FeatureCurator().curate(entry, feature));
    }

    @Test
    public void canKeepTrackOfFeatures() {
        FeatureCurator curator = new FeatureCurator();

        ViccEntry entry = ViccTestFactory.testViccEntryForSource(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript("any"));
        Feature feature = ImmutableFeature.builder().geneSymbol("any").name("any").build();

        assertNotNull(curator.curate(entry, feature));
        int unusedCurationKeyCount = curator.unusedCurationKeys().size();

        CurationKey blacklistKey = CurationFactory.ONCOKB_FEATURE_BLACKLIST.iterator().next();
        ViccEntry blacklistEntry =
                ViccTestFactory.testViccEntryForSource(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript(blacklistKey.transcript()));

        Feature blacklistFeature = ImmutableFeature.builder().geneSymbol(blacklistKey.gene()).name(blacklistKey.featureName()).build();

        assertNull(curator.curate(blacklistEntry, blacklistFeature));
        int newUnusedCurationKeyCount = curator.unusedCurationKeys().size();
        assertEquals(1, unusedCurationKeyCount - newUnusedCurationKeyCount);
    }
}