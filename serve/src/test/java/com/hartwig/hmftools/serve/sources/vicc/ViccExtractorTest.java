package com.hartwig.hmftools.serve.sources.vicc;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.genome.genepanel.HmfGenePanelSupplier;
import com.hartwig.hmftools.common.genome.region.HmfTranscriptRegion;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.extraction.hotspot.ProteinResolver;
import com.hartwig.hmftools.vicc.annotation.ViccClassificationConfig;
import com.hartwig.hmftools.vicc.datamodel.Association;
import com.hartwig.hmftools.vicc.datamodel.ViccEntry;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class ViccExtractorTest {

    private static final Map<String, HmfTranscriptRegion> V37_GENE_MAP = HmfGenePanelSupplier.allGenesMap37();

    @Test
    public void canExtractFromViccEntries() throws IOException {
        EventClassifierConfig config = ViccClassificationConfig.build();
        ViccExtractor extractor =
                ViccExtractorFactory.buildViccExtractor(config, new TestProteinResolver(), Lists.newArrayList(), V37_GENE_MAP);

        Association association =
                ViccTestFactory.testActionableAssociation("drugs", "colorectal cancer", "DOID:123", "A", "Responsive", "http");

        List<ViccEntry> entries = Lists.newArrayList();
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("KIT", "KIT Amplification", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("BRAF", "V600E", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("NTRK3", "NTRK3 Fusion", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("BRAF", "V600", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("EGFR", "Exon 19 deletion", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("ALK", "EML4-ALK Fusion", association));
        entries.add(ViccTestFactory.testEntryWithGeneEventAndAssociation("-", "Microsatellite Instability-High", association));

        ExtractionResult result = extractor.extract(entries);
        assertEquals(1, result.knownHotspots().size());
        assertEquals(1, result.knownCopyNumbers().size());
        assertEquals(1, result.knownFusionPairs().size());
        assertEquals(1, result.actionableHotspots().size());
        assertEquals(2, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(1, result.actionableSignatures().size());
    }

    private static class TestProteinResolver implements ProteinResolver {

        @NotNull
        @Override
        public List<VariantHotspot> resolve(@NotNull final String gene, @Nullable final String specificTranscript,
                @NotNull final String proteinAnnotation) {
            return Lists.newArrayList(ImmutableVariantHotspotImpl.builder().chromosome("1").position(10).ref("A").alt("T").build());
        }

        @NotNull
        @Override
        public Set<String> unresolvedProteinAnnotations() {
            return Sets.newHashSet();
        }
    }
}