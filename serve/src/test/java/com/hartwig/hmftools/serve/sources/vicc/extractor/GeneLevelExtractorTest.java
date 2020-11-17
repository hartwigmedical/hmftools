package com.hartwig.hmftools.serve.sources.vicc.extractor;

import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.ONCO;
import static com.hartwig.hmftools.common.drivercatalog.DriverCategory.TSG;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.serve.actionability.gene.GeneLevelEvent;
import com.hartwig.hmftools.vicc.datamodel.Feature;
import com.hartwig.hmftools.vicc.datamodel.ImmutableFeature;

import org.jetbrains.annotations.NotNull;
import org.junit.Ignore;
import org.junit.Test;

public class GeneLevelExtractorTest {

    @Test
    @Ignore
    public void canExtractGeneLevelEvent() {
        //TODO improve test
        Feature feature = ImmutableFeature.builder().geneSymbol("a").name("a").build();

        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(feature, createDriverGeneONCO()), GeneLevelEvent.ACTIVATION);
        assertEquals(GeneLevelExtractor.extractGeneLevelEvent(feature, createDriverGeneTSG()), GeneLevelEvent.INACTIVATION);
    }

    @NotNull
    private static List<DriverGene> createDriverGeneTSG() {
        return Lists.newArrayList(ImmutableDriverGene.builder()
                .gene("a")
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(TSG)
                .build());
    }

    @NotNull
    private static List<DriverGene> createDriverGeneONCO() {
        return Lists.newArrayList(ImmutableDriverGene.builder()
                .gene("a")
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(true)
                .reportSomaticHotspot(true)
                .reportGermlineVariant(false)
                .reportGermlineHotspot(false)
                .likelihoodType(ONCO)
                .build());
    }
}