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

public class ViccCuratorTest {

    @Test
    public void canCurateFeatures() {
        CurationKey firstOncoKbKey = firstOncoKbMappingKey();
        String firstMappedFeature = CurationFactory.FEATURE_MAPPINGS.get(firstOncoKbKey).featureName();

        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                firstOncoKbKey.transcript(),
                "gene",
                "event",
                "description",
                "chromosome",
                "pos",
                null);

        Feature feature = ViccTestFactory.testFeatureWithGeneAndName(firstOncoKbKey.gene(),
                firstOncoKbKey.featureName(),
                "description",
                "chromosome",
                "pos",
                null);

        assertEquals(firstMappedFeature, new ViccCurator().curate(entry, feature).name());
    }

    @Test
    public void canBlacklistFeatures() {
        CurationKey firstOncoKbKey = firstOncoKbBlacklistKey();
        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                firstOncoKbKey.transcript(),
                "gene",
                "event",
                "description",
                "chromosome",
                "pos",
                null);

        Feature feature = ViccTestFactory.testFeatureWithGeneAndName(firstOncoKbKey.gene(),
                firstOncoKbKey.featureName(),
                "description",
                "chromosome",
                "pos",
                null);
        assertNull(new ViccCurator().curate(entry, feature));
    }

    @Test
    public void canKeepTrackOfFeatures() {
        ViccCurator curator = new ViccCurator();

        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                "any",
                "gene",
                "event",
                "description",
                "chromosome",
                "pos",
                null);
        Feature feature = ViccTestFactory.testFeatureWithGeneAndName("any", "description", "any", "chromosome", "pos", null);

        assertNotNull(curator.curate(entry, feature));

        CurationKey blacklistKey = firstOncoKbBlacklistKey();
        ViccEntry blacklistEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                blacklistKey.transcript(),
                "gene",
                "event",
                "description",
                "chromosome",
                "pos",
                null);

        Feature blacklistFeature = ViccTestFactory.testFeatureWithGeneAndName(blacklistKey.gene(),
                blacklistKey.featureName(),
                "description",
                "chromosome",
                "pos",
                null);

        assertNull(curator.curate(blacklistEntry, blacklistFeature));

        curator.reportUnusedCurationEntries();
    }

    @NotNull
    private static CurationKey firstOncoKbMappingKey() {
        for (CurationKey key : CurationFactory.FEATURE_MAPPINGS.keySet()) {
            if (key.source() == ViccSource.ONCOKB) {
                return key;
            }
        }
        throw new IllegalStateException("No OncoKB mapping keys found!");
    }

    @NotNull
    private static CurationKey firstOncoKbBlacklistKey() {
        for (CurationKey key : CurationFactory.FEATURE_BLACKLIST) {
            if (key.source() == ViccSource.ONCOKB) {
                return key;
            }
        }
        throw new IllegalStateException("No OncoKB blacklist keys found!");
    }
}