package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.serve.sources.vicc.ViccTestFactory;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class HotspotExtractorTest {

    private static final VariantHotspot TEST_HOTSPOT =
            ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();

    @Test
    public void canExtractHotspot() {
        Feature hotspotFeature = ViccTestFactory.testFeatureWithName("V600E");
        Feature ampFeature = ViccTestFactory.testFeatureWithName("Amplification");
        ViccEntry entry = ViccTestFactory.testEntryWithFeatures(Lists.newArrayList(hotspotFeature, ampFeature));

        HotspotExtractor hotspotExtractor = new HotspotExtractor(new TestProteinResolver(), new ProteinAnnotationExtractor());
        Map<Feature, List<VariantHotspot>> hotspots = hotspotExtractor.extractHotspots(entry);

        assertEquals(1, hotspots.size());
        assertEquals(1, hotspots.get(hotspotFeature).size());
        assertEquals(TEST_HOTSPOT, hotspots.get(hotspotFeature).get(0));
        assertNull(hotspots.get(ampFeature));
    }

    private static class TestProteinResolver implements ProteinResolver {

        @NotNull
        @Override
        public List<VariantHotspot> extractHotspotsFromProteinAnnotation(@NotNull final String gene,
                @Nullable final String specificTranscript, @NotNull final String proteinAnnotation) {
            return Lists.newArrayList(TEST_HOTSPOT);
        }

        @NotNull
        @Override
        public Set<String> unresolvedProteinAnnotations() {
            return Sets.newHashSet();
        }
    }
}