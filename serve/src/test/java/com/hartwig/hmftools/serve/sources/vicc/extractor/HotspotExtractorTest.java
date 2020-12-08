package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.GeneChecker;
import com.hartwig.hmftools.serve.extraction.GeneCheckerTestFactory;
import com.hartwig.hmftools.serve.hotspot.ProteinResolver;
import com.hartwig.hmftools.vicc.annotation.ProteinAnnotationExtractor;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class HotspotExtractorTest {

    private static final GeneChecker HG19_GENE_CHECKER = GeneCheckerTestFactory.buildForHG19();

    private static final VariantHotspot TEST_HOTSPOT =
            ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();

    @Test
    public void canExtractSimpleHotspot() {
        String protein = "V600E";

        HotspotExtractor hotspotExtractor = createWithProtein(protein);
        List<VariantHotspot> hotspots = hotspotExtractor.extract("BRAF", null, EventType.HOTSPOT, "V600E");

        assertEquals(1, hotspots.size());
        assertEquals(TEST_HOTSPOT, hotspots.get(0));
    }

    @Test
    public void skipsHotspotsForInvalidGenes() {
        String protein = "V600E";

        HotspotExtractor hotspotExtractor = createWithProtein(protein);

        assertNull(hotspotExtractor.extract("NOT-A-GENE", null, EventType.HOTSPOT, "V600E"));
    }

    @NotNull
    private static HotspotExtractor createWithProtein(@NotNull String protein) {
        return new HotspotExtractor(HG19_GENE_CHECKER, new TestProteinResolver(protein), new ProteinAnnotationExtractor());
    }

    private static class TestProteinResolver implements ProteinResolver {

        @NotNull
        private final String protein;

        public TestProteinResolver(@NotNull final String protein) {
            this.protein = protein;
        }

        @NotNull
        @Override
        public List<VariantHotspot> resolve(@NotNull final String gene, @Nullable final String specificTranscript,
                @NotNull final String proteinAnnotation) {
            return proteinAnnotation.equals(protein) ? Lists.newArrayList(TEST_HOTSPOT) : Lists.newArrayList();
        }

        @NotNull
        @Override
        public Set<String> unresolvedProteinAnnotations() {
            return Sets.newHashSet();
        }
    }
}