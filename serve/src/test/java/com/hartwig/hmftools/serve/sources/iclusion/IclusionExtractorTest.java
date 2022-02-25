package com.hartwig.hmftools.serve.sources.iclusion;

import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.or;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTumorLocation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutation;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionTumorLocation;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;

import org.junit.Test;

public class IclusionExtractorTest {

    @Test
    public void canExtractFromIclusionEntries() {
        EventClassifierConfig config = IclusionClassificationConfig.build();
        IclusionExtractor extractor = IclusionExtractorFactory.buildIclusionExtractor(config,
                RefGenomeResourceTestFactory.buildTestResource37(),
                DoidLookupTestFactory.dummy());

        IclusionTumorLocation loc1 =
                ImmutableIclusionTumorLocation.builder().primaryTumorLocation("ptum").addDoids("162").build();

        List<IclusionTrial> entries = Lists.newArrayList();
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("KIT")
                        .name("KIT AMPLIFICATION")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("BRAF")
                        .name("V600E")
                        .negation(false)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("NTRK3")
                        .name("NTRK3 FUSION")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("BRAF")
                        .name("V600")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("BRAF")
                        .name("EXON 1 DELETION")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("ALK")
                        .name("EML4-ALK Fusion")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));
        entries.add(IclusionTestFactory.trialWithMutationsAndTumorLocation("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("-")
                        .name("MSI HIGH")
                        .negation(true)
                        .build()))), Lists.newArrayList(loc1)));

        ExtractionResult result = extractor.extract(entries);
        // Iclusion don't extract known events
        assertEquals(0, result.knownHotspots().size());
        assertEquals(0, result.knownCopyNumbers().size());
        assertEquals(0, result.knownFusionPairs().size());

        assertEquals(1, result.actionableHotspots().size());
        assertEquals(2, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(1, result.actionableCharacteristics().size());
    }
}