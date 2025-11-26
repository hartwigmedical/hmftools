package com.hartwig.hmftools.orange.algo.linx;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.datamodel.finding.Fusion;

import org.junit.Test;

public class DnaFusionEvaluatorTest
{
    @Test
    public void canDetermineIfFusionIsPresent()
    {
        Fusion fusion1 = LinxOrangeTestFactory.fusionBuilder().geneStart("start 1").geneEnd("end 1").build();
        Fusion fusion2 = LinxOrangeTestFactory.fusionBuilder().geneStart("start 2").geneEnd("end 1").build();

        List<Fusion> fusions = Lists.newArrayList(fusion1, fusion2);
        assertTrue(DnaFusionEvaluator.hasFusion(fusions, "start 2", "end 1"));
        assertFalse(DnaFusionEvaluator.hasFusion(fusions, "start 2", "end 2"));
        assertFalse(DnaFusionEvaluator.hasFusion(fusions, "start 1", "end 2"));
    }
}