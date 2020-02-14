package com.hartwig.hmftools.svtools.rna_expression;

import static com.hartwig.hmftools.linx.types.SvVarData.SE_END;
import static com.hartwig.hmftools.linx.types.SvVarData.SE_START;
import static com.hartwig.hmftools.svtools.rna_expression.ExpectedExpressionRates.readsSupportFragment;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.variant.structural.annotation.ExonData;
import com.hartwig.hmftools.common.variant.structural.annotation.TranscriptData;

import org.junit.Test;

public class RnaExpectedRates
{
    @Test
    public void testRegionMatching()
    {
        RnaExpConfig config = new RnaExpConfig();
        config.MedianFragmentLength = 100;
        config.ReadLength = 20;

        ExpectedExpressionRates eeRates = new ExpectedExpressionRates(config);

        TranscriptData transData = new TranscriptData(1, "TRANS01", "GENE01", true, (byte)1,
                0, 1000, null,null, "");

        transData.exons().add(new ExonData(1, 100, 200, 1, -1, -1));
        transData.exons().add(new ExonData(1, 300, 400, 2, -1, -1));
        transData.exons().add(new ExonData(1, 440, 449, 3, -1, -1));
        transData.exons().add(new ExonData(1, 460, 469, 4, -1, -1));
        transData.exons().add(new ExonData(1, 600, 800, 5, -1, -1));

        // fully contained fragment
        long startPos = 100;
        List<long[]> readRegions = Lists.newArrayList();
        boolean spansExons = eeRates.generateImpliedFragment(transData, startPos, readRegions);

        assertFalse(spansExons);
        assertEquals(2, readRegions.size());
        assertEquals(100, readRegions.get(0)[SE_START]);
        assertEquals(119, readRegions.get(0)[SE_END]);
        assertEquals(180, readRegions.get(1)[SE_START]);
        assertEquals(199, readRegions.get(1)[SE_END]);
        assertTrue(readsSupportFragment(transData, readRegions, spansExons));

        // 2 fully contained reads but in 2 exons
        startPos = 371;
        spansExons = eeRates.generateImpliedFragment(transData, startPos, readRegions);

        assertTrue(spansExons);
        assertEquals(2, readRegions.size());
        assertEquals(371, readRegions.get(0)[SE_START]);
        assertEquals(390, readRegions.get(0)[SE_END]);
        assertEquals(631, readRegions.get(1)[SE_START]);
        assertEquals(650, readRegions.get(1)[SE_END]);
        assertTrue(readsSupportFragment(transData, readRegions, spansExons));

        // within an exon then spanning a junction

        startPos = 150;
        config.MedianFragmentLength = 32 + 85 + 2 * config.ReadLength;
        spansExons = eeRates.generateImpliedFragment(transData, startPos, readRegions);

        assertTrue(spansExons);
        assertEquals(3, readRegions.size());
        assertEquals(150, readRegions.get(0)[SE_START]);
        assertEquals(169, readRegions.get(0)[SE_END]);
        assertEquals(387, readRegions.get(1)[SE_START]);
        assertEquals(400, readRegions.get(1)[SE_END]);
        assertEquals(440, readRegions.get(2)[SE_START]);
        assertEquals(445, readRegions.get(2)[SE_END]);
        assertTrue(readsSupportFragment(transData, readRegions, spansExons));

        // start spanning a junction, ends within an exon
        startPos = 191;
        config.MedianFragmentLength = 81 + 20 + 15 + 2 * config.ReadLength;
        spansExons = eeRates.generateImpliedFragment(transData, startPos, readRegions);

        assertTrue(spansExons);
        assertEquals(3, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(606, readRegions.get(2)[SE_START]);
        assertEquals(625, readRegions.get(2)[SE_END]);
        assertTrue(readsSupportFragment(transData, readRegions, spansExons));

        // 2 sets of exon junctions
        startPos = 191;
        config.MedianFragmentLength = 81 + 5 + 2 * config.ReadLength;
        spansExons = eeRates.generateImpliedFragment(transData, startPos, readRegions);

        assertTrue(spansExons);
        assertEquals(5, readRegions.size());
        assertEquals(191, readRegions.get(0)[SE_START]);
        assertEquals(200, readRegions.get(0)[SE_END]);
        assertEquals(300, readRegions.get(1)[SE_START]);
        assertEquals(309, readRegions.get(1)[SE_END]);
        assertEquals(396, readRegions.get(2)[SE_START]);
        assertEquals(400, readRegions.get(2)[SE_END]);
        assertEquals(440, readRegions.get(3)[SE_START]);
        assertEquals(449, readRegions.get(3)[SE_END]);
        assertEquals(460, readRegions.get(4)[SE_START]);
        assertEquals(464, readRegions.get(4)[SE_END]);
        assertTrue(readsSupportFragment(transData, readRegions, spansExons));
    }

}
