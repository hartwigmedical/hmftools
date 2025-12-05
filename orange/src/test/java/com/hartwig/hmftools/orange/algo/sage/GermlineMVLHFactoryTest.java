package com.hartwig.hmftools.orange.algo.sage;

import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO;
import static com.hartwig.hmftools.common.purple.PurpleCommon.DEFAULT_DRIVER_HET_DELETION_THRESHOLD;

import static org.apache.commons.math3.util.Precision.EPSILON;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotNull;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.driver.DriverCategory;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.driver.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.driver.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.metrics.GeneDepth;

import org.junit.Test;

public class GermlineMVLHFactoryTest
{
    @Test
    public void parseGermlineMVLHCorrectly()
    {
        String reliableMVLHGene = "FAKE1";
        String unReliableMVLHGene = "FAKE2";

        final int[] emptyDepth = { 0, 0 };

        List<GeneDepth> geneDepths = List.of(
                new GeneDepth(reliableMVLHGene, "7", 1000, 2000, 0.0002, emptyDepth),
                new GeneDepth(unReliableMVLHGene, "1", 1000, 2000, 1.0000, emptyDepth));

        DriverGene driverGene1 = ImmutableDriverGene.builder()
                .gene(reliableMVLHGene)
                .reportMissenseAndInframe(true)
                .reportNonsenseAndFrameshift(true)
                .reportSplice(true)
                .reportDeletion(true)
                .reportDisruption(true)
                .reportAmplification(false)
                .amplificationRatio(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                .reportHetDeletion(false)
                .hetDeletionThreshold(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                .reportSomaticHotspot(true)
                .likelihoodType(DriverCategory.TSG)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build();
        DriverGene driverGene2 = ImmutableDriverGene.builder()
                .gene(unReliableMVLHGene)
                .reportMissenseAndInframe(false)
                .reportNonsenseAndFrameshift(false)
                .reportSplice(false)
                .reportDeletion(false)
                .reportDisruption(false)
                .reportAmplification(true)
                .amplificationRatio(DEFAULT_DRIVER_AMPLIFICATION_PLOIDY_RATIO)
                .reportHetDeletion(false)
                .hetDeletionThreshold(DEFAULT_DRIVER_HET_DELETION_THRESHOLD)
                .reportSomaticHotspot(false)
                .likelihoodType(DriverCategory.ONCO)
                .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                .reportGermlineDeletion(DriverGeneGermlineReporting.NONE)
                .reportGermlineDisruption(DriverGeneGermlineReporting.NONE)
                .reportPGX(false)
                .build();

        Map<String,DriverGene> driverGenes = Maps.newHashMap();
        driverGenes.put(driverGene1.gene(), driverGene1);
        driverGenes.put(driverGene2.gene(), driverGene2);

        Map<String, Double> mvlhPerGene = GermlineMVLHFactory.parseMVLHPerGene(geneDepths, driverGenes);

        assertNotNull(mvlhPerGene);
        assertEquals(0.0002D, mvlhPerGene.get(reliableMVLHGene), EPSILON);
        assertFalse(mvlhPerGene.containsKey(unReliableMVLHGene));
    }
}