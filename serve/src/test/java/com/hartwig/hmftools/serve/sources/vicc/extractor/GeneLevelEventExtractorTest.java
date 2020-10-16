package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.*;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.junit.Ignore;
import org.junit.Test;

public class GeneLevelEventExtractorTest {

//    @Ignore
//    @Test
//    public void canExtractGeneLevelEvent() {
//        Feature feature = ImmutableFeature.builder().geneSymbol("a").name("a").build();
//
//        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("promiscuousFusion", feature, createDriverGene("Gene1")), GeneLevelEvent.FUSION);
//        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("geneLevel", feature, createDriverGene("Gene1")), GeneLevelEvent.ACTIVATION);
//        assertEquals(GeneLevelEventExtractor.extractGeneLevelEvent("ccc", feature, createDriverGene("Gene1")), GeneLevelEvent.UNKONWN);
//        assertNotEquals(GeneLevelEventExtractor.extractGeneLevelEvent("aa", feature, createDriverGene("Gene1")), GeneLevelEvent.INACTIVATION);
//    }

    private static List<DriverGene> createDriverGene(final String name)
    {
        return Lists.newArrayList(ImmutableDriverGene.builder()
                .gene(name)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(true)
                .reportAmplification(false)
                .reportHotspot(false)
                .likelihoodType(TSG)
                .build());
    }

}