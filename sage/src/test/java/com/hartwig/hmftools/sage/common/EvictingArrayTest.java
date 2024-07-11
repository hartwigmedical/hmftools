package com.hartwig.hmftools.sage.common;

import static com.hartwig.hmftools.common.hla.HlaCommon.HLA_GENES;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.GENE_ID_1;
import static com.hartwig.hmftools.sage.SageConstants.REGION_BLOCK_SIZE;
import static com.hartwig.hmftools.sage.common.TestUtils.region;
import static com.hartwig.hmftools.sage.select.ReadPanelStatus.MIXED;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.function.Consumer;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.genome.region.Strand;
import com.hartwig.hmftools.common.hla.HlaCommon;
import com.hartwig.hmftools.common.region.BasePosition;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.sage.candidate.EvictingArray;
import com.hartwig.hmftools.sage.candidate.RefContext;
import com.hartwig.hmftools.sage.candidate.RegionBlock;
import com.hartwig.hmftools.sage.select.PanelSelector;

import org.junit.Test;

public class EvictingArrayTest
{
    @Test
    public void testRegionBlocks()
    {
        List<BaseRegion> panel = Lists.newArrayList(region(950, 1050), region(1100, 1450));

        PanelSelector panelSelector = new PanelSelector(panel);

        String chromosome = "6";

        HlaCommon.populateGeneData(List.of(new GeneData(
                GENE_ID_1, HLA_GENES.get(0), chromosome, Strand.POS_STRAND, 1000, 1300, "")));

        ChrBaseRegion chrBaseRegion = new ChrBaseRegion(chromosome, 700, 1500);

        List<SimpleVariant> hotspots = Lists.newArrayList(
                new SimpleVariant(chromosome, 910, "A", "T"),
                new SimpleVariant(chromosome, 1350, "A", "T"));

        List<RegionBlock> regionBlocks = RegionBlock.buildRegionBlocks(
                REGION_BLOCK_SIZE, chrBaseRegion, panelSelector, hotspots, 20, 10);

        assertEquals(8, regionBlocks.size());

        assertEquals(10, regionBlocks.get(0).depthLimit());
        assertEquals(10, regionBlocks.get(1).depthLimit());
        assertEquals(20, regionBlocks.get(2).depthLimit());
        assertEquals(MIXED, regionBlocks.get(2).panelStatus());
        assertEquals(20, regionBlocks.get(7).depthLimit());

        assertTrue(regionBlocks.get(0).applyEventPenalty());
        assertFalse(regionBlocks.get(4).applyEventPenalty());

        assertTrue(regionBlocks.get(2).coversHotspot(900, 920));

        RegionBlock regionBlock = regionBlocks.get(0);

        for(int i = 0; i < 9; ++i)
        {
            regionBlock.incrementDepth();
        }

        assertFalse(regionBlock.depthLimitReached());
        regionBlock.incrementDepth();
        assertTrue(regionBlock.depthLimitReached());
    }

    @Test
    public void testEvictingArray()
    {
        int readLength = 200;

        EvictionHandler handler = new EvictionHandler();
        EvictingArray evictingArray = new EvictingArray(readLength, handler);

        int capacity = evictingArray.capacity();
        assertEquals(readLength * 2, capacity);

        int startPosition = 1000;
        int minPosition = startPosition - evictingArray.readLengthBuffer();
        evictingArray.getOrCreateRefContext(startPosition, EvictingArrayTest::create);
        assertEquals(minPosition, evictingArray.minPosition());

        int position = startPosition;
        for(; position < minPosition + capacity; position++)
        {
            evictingArray.getOrCreateRefContext(position, EvictingArrayTest::create);
        }

        assertEquals(120, evictingArray.itemCount());

        // no flush yet
        assertEquals(0, handler.items().size());

        for(; position < 1400; position++)
        {
            evictingArray.getOrCreateRefContext(position, EvictingArrayTest::create);
        }

        assertEquals(1080, evictingArray.minPosition());
        assertEquals(320, evictingArray.itemCount());
        assertEquals(80, handler.items().size());

        for(; position < 1600; position++)
        {
            evictingArray.getOrCreateRefContext(position, EvictingArrayTest::create);
        }

        assertEquals(1200, evictingArray.minPosition());
        assertEquals(400, evictingArray.itemCount());
        assertEquals(200, handler.items().size());

        // evict all but the last

        evictingArray.getOrCreateRefContext(2000, EvictingArrayTest::create);
        assertEquals(1, evictingArray.itemCount());
        assertEquals(600, handler.items().size());


        evictingArray.evictAll();
        assertEquals(0, evictingArray.itemCount());
        assertEquals(601, handler.items().size());
    }

    static class EvictionHandler implements Consumer<RefContext>
    {
        private final List<BasePosition> mItems = Lists.newArrayList();

        public List<BasePosition> items() { return mItems; }

        @Override
        public void accept(final RefContext position)
        {
            mItems.add(position);
        }
    }

    private static RefContext create(int pos)
    {
        return new RefContext(CHR_1, pos);
    }
}
