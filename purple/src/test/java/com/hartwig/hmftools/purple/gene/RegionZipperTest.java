package com.hartwig.hmftools.purple.gene;

import static org.mockito.Mockito.mock;
import static org.mockito.Mockito.verify;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.region.GenomeRegion;
import com.hartwig.hmftools.common.genome.region.GenomeRegions;

import org.junit.Test;

public class RegionZipperTest {

    @Test
    public void testZipper() {
        @SuppressWarnings("unchecked")
        RegionZipperHandler<GenomeRegion, GenomeRegion> handler = mock(RegionZipperHandler.class);

        final List<GenomeRegion> primary = Lists.newArrayList(GenomeRegions.create("1", 1, 100),
                GenomeRegions.create("1", 500, 1000),
                GenomeRegions.create("2", 500, 1000),
                GenomeRegions.create("4", 500, 1000));

        final List<GenomeRegion> secondary = Lists.newArrayList(
                GenomeRegions.create("1", 500, 1000),
                GenomeRegions.create("2", 500, 1000),
                GenomeRegions.create("3", 500, 1000));

        RegionZipper.zip(primary, secondary, handler);

        verify(handler).enterChromosome("1");
        verify(handler).primary(primary.get(0));
        verify(handler).primary(primary.get(1));
        verify(handler).secondary(secondary.get(0));
        verify(handler).enterChromosome("2");
        verify(handler).primary(primary.get(2));
        verify(handler).secondary(secondary.get(1));
        verify(handler).enterChromosome("3");
        verify(handler).secondary(secondary.get(2));
        verify(handler).enterChromosome("4");
        verify(handler).primary(primary.get(3));
    }
}
