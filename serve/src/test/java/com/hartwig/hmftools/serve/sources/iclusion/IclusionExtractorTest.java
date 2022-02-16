package com.hartwig.hmftools.serve.sources.iclusion;

import static com.hartwig.hmftools.serve.sources.iclusion.IclusionTestFactory.or;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.serve.classification.EventClassifierConfig;
import com.hartwig.hmftools.iclusion.classification.IclusionClassificationConfig;
import com.hartwig.hmftools.iclusion.datamodel.IclusionTrial;
import com.hartwig.hmftools.iclusion.datamodel.ImmutableIclusionMutation;
import com.hartwig.hmftools.serve.curation.DoidLookupTestFactory;
import com.hartwig.hmftools.serve.extraction.ExtractionResult;
import com.hartwig.hmftools.serve.refgenome.RefGenomeResourceTestFactory;

import org.junit.Ignore;
import org.junit.Test;

public class IclusionExtractorTest {

    @Test
    @Ignore
    public void canExtractFromIclusionEntries() {
        EventClassifierConfig config = IclusionClassificationConfig.build();
        IclusionExtractor extractor = IclusionExtractorFactory.buildIclusionExtractor(config,
                RefGenomeResourceTestFactory.buildTestResource37(),
                DoidLookupTestFactory.dummy());

        List<IclusionTrial> entries = Lists.newArrayList();
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("KIT")
                        .name("KIT AMPLIFICATION")
                        .negation(true)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("BRAF")
                        .name("V600E")
                        .negation(false)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("NTRK3")
                        .name("NTRK3 FUSION")
                        .negation(true)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("BRAF")
                        .name("V600")
                        .negation(true)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("EGFR")
                        .name("Exon 19")
                        .negation(true)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("ALK")
                        .name("EML4-ALK Fusion")
                        .negation(true)
                        .build())))));
        entries.add(IclusionTestFactory.trial("trial",
                Lists.newArrayList(or(Lists.newArrayList(ImmutableIclusionMutation.builder()
                        .gene("-")
                        .name("MSI HIGH")
                        .negation(true)
                        .build())))));

        ExtractionResult result = extractor.extract(entries);
        assertEquals(1, result.knownHotspots().size());
        assertEquals(1, result.knownCopyNumbers().size());
        assertEquals(1, result.knownFusionPairs().size());
        assertEquals(1, result.actionableHotspots().size());
        assertEquals(2, result.actionableRanges().size());
        assertEquals(2, result.actionableGenes().size());
        assertEquals(1, result.actionableFusions().size());
        assertEquals(1, result.actionableCharacteristics().size());

    }
}