package com.hartwig.hmftools.serve.extraction.hotspot;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.serve.classification.EventType;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.DriverGenesTestFactory;
import com.hartwig.hmftools.serve.extraction.util.DriverInconsistencyMode;
import com.hartwig.hmftools.serve.extraction.util.GeneChecker;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class HotspotExtractorTest {

    private static final GeneChecker GENE_CHECKER = new GeneChecker(Sets.newHashSet("BRAF", "KRAS"));

    private static final VariantHotspot TEST_HOTSPOT =
            ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build();

    @Test
    public void canFindGene() {
        List<DriverGene> driverGenes = DriverGenesTestFactory.createDriverGenes("BRAF", "KIT");
        assertEquals(DriverCategory.TSG, HotspotExtractor.findByGene(driverGenes, "BRAF"));
        assertEquals(DriverCategory.ONCO, HotspotExtractor.findByGene(driverGenes, "KIT"));
        assertNull(HotspotExtractor.findByGene(driverGenes, "KRAS"));
    }

    @Test
    public void canFilterInCatalog() {
        String protein = "V600E";
        HotspotExtractor hotspotExtractorFilter = createWithProtein(protein, DriverInconsistencyMode.FILTER);
        List<VariantHotspot> hotspotExtractorFilterList = hotspotExtractorFilter.extract("BRAF", null, EventType.HOTSPOT, "V600E");
        assertEquals(1, hotspotExtractorFilterList.size());
        assertEquals(TEST_HOTSPOT, hotspotExtractorFilterList.get(0));

        HotspotExtractor hotspotExtractorIgnore = createWithProtein(protein, DriverInconsistencyMode.IGNORE);
        List<VariantHotspot> hotspotExtractorIgnoreList = hotspotExtractorIgnore.extract("BRAF", null, EventType.HOTSPOT, "V600E");
        assertEquals(1, hotspotExtractorIgnoreList.size());
        assertEquals(TEST_HOTSPOT, hotspotExtractorIgnoreList.get(0));

        HotspotExtractor hotspotExtractorWarn = createWithProtein(protein, DriverInconsistencyMode.WARN_ONLY);
        List<VariantHotspot> hotspotExtractorWarnList = hotspotExtractorWarn.extract("BRAF", null, EventType.HOTSPOT, "V600E");
        assertEquals(1, hotspotExtractorWarnList.size());
        assertEquals(TEST_HOTSPOT, hotspotExtractorWarnList.get(0));
    }

    @Test
    public void canFilterNotInCatalog() {
        String protein = "V600E";
        HotspotExtractor hotspotExtractorFilter = createWithProtein(protein, DriverInconsistencyMode.FILTER);
        List<VariantHotspot> hotspotExtractorFilterList = hotspotExtractorFilter.extract("KRAS", null, EventType.HOTSPOT, "V600E");
        assertNull(hotspotExtractorFilterList);

        HotspotExtractor hotspotExtractorIgnore = createWithProtein(protein, DriverInconsistencyMode.IGNORE);
        List<VariantHotspot> hotspotExtractorIgnoreList = hotspotExtractorIgnore.extract("KRAS", null, EventType.HOTSPOT, "V600E");
        assertEquals(1, hotspotExtractorIgnoreList.size());
        assertEquals(TEST_HOTSPOT, hotspotExtractorIgnoreList.get(0));

        HotspotExtractor hotspotExtractorWarn = createWithProtein(protein, DriverInconsistencyMode.WARN_ONLY);
        List<VariantHotspot> hotspotExtractorWarnList = hotspotExtractorWarn.extract("KRAS", null, EventType.HOTSPOT, "V600E");
        assertEquals(1, hotspotExtractorWarnList.size());
        assertEquals(TEST_HOTSPOT, hotspotExtractorWarnList.get(0));
    }

    @Test
    public void canExtractSimpleHotspot() {
        String protein = "V600E";

        HotspotExtractor hotspotExtractor = createWithProtein(protein, DriverInconsistencyMode.IGNORE);
        List<VariantHotspot> hotspots = hotspotExtractor.extract("BRAF", null, EventType.HOTSPOT, "V600E");

        assertEquals(1, hotspots.size());
        assertEquals(TEST_HOTSPOT, hotspots.get(0));
    }

    @Test
    public void skipsHotspotsForInvalidGenes() {
        String protein = "V600E";

        HotspotExtractor hotspotExtractor = createWithProtein(protein, DriverInconsistencyMode.IGNORE);

        assertNull(hotspotExtractor.extract("NOT-A-GENE", null, EventType.HOTSPOT, "V600E"));
    }

    @NotNull
    private static HotspotExtractor createWithProtein(@NotNull String protein, @NotNull DriverInconsistencyMode annotation) {
        return new HotspotExtractor(GENE_CHECKER,
                new TestProteinResolver(protein),
                event -> event,
                annotation,
                DriverGenesTestFactory.createDriverGenes("BRAF", "KIT"));
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