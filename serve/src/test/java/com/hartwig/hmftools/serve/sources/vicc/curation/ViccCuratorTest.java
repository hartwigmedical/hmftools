package com.hartwig.hmftools.serve.sources.vicc.curation;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;
import com.hartwig.hmftools.vicc.datamodel.ViccSource;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class ViccCuratorTest {

    @Test
    public void canCurateFeatures() {
        CurationKey firstOncoKbKey = firstOncoKbMappingKey();
        String firstMappedFeature = CurationFactory.FEATURE_NAME_MAPPINGS.get(firstOncoKbKey);

        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                ViccTestFactory.testOncoKbWithTranscript(firstOncoKbKey.transcript()));

        Feature feature = ImmutableFeature.builder().geneSymbol(firstOncoKbKey.gene()).name(firstOncoKbKey.featureName()).build();

        assertEquals(firstMappedFeature, new ViccCurator().curate(entry, feature).name());
    }

    @Test
    public void canBlacklistFeatures() {
        CurationKey firstOncoKbKey = firstOncoKbBlacklistKey();
        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                ViccTestFactory.testOncoKbWithTranscript(firstOncoKbKey.transcript()));

        Feature feature = ImmutableFeature.builder().geneSymbol(firstOncoKbKey.gene()).name(firstOncoKbKey.featureName()).build();
        assertNull(new ViccCurator().curate(entry, feature));
    }

    @Test
    @Ignore
    public void canKeepTrackOfFeatures() {
        ViccCurator curator = new ViccCurator();

        ViccEntry entry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB, ViccTestFactory.testOncoKbWithTranscript("any"));
        Feature feature = ImmutableFeature.builder().geneSymbol("any").name("any").build();

        assertNotNull(curator.curate(entry, feature));

        CurationKey blacklistKey = firstOncoKbBlacklistKey();
        ViccEntry blacklistEntry = ViccTestFactory.testViccEntryWithSourceAndKbObject(ViccSource.ONCOKB,
                ViccTestFactory.testOncoKbWithTranscript(blacklistKey.transcript()));

        Feature blacklistFeature = ImmutableFeature.builder().geneSymbol(blacklistKey.gene()).name(blacklistKey.featureName()).build();

        assertNull(curator.curate(blacklistEntry, blacklistFeature));

        curator.reportUnusedCurationEntries();
    }

    @NotNull
    private static CurationKey firstOncoKbMappingKey() {
        for (CurationKey key : CurationFactory.FEATURE_NAME_MAPPINGS.keySet()) {
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