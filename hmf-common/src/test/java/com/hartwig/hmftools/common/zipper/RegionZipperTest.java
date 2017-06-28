package com.hartwig.hmftools.common.zipper;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.bed.ImmutableBEDGenomeRegion;

import org.junit.Test;

import mockit.Expectations;
import mockit.Injectable;

public class RegionZipperTest {

    @Injectable
    private RegionZipperHandler<GenomeRegion, GenomeRegion> handler;

    @Test
    public void testZipper() {

        final List<GenomeRegion> primary = Lists.newArrayList(
                createRegion("1", 1, 100),
                createRegion("1", 500, 1000),
                createRegion("2", 500, 1000),
                createRegion("4", 500, 1000));

        final List<GenomeRegion> secondary = Lists.newArrayList(
                createRegion("1", 500, 1000),
                createRegion("2", 500, 1000),
                createRegion("3", 500, 1000));

        new Expectations() {{
            handler.enterChromosome("1");
            handler.primary(primary.get(0));
            handler.primary(primary.get(1));
            handler.secondary(secondary.get(0));
            handler.enterChromosome("2");
            handler.primary(primary.get(2));
            handler.secondary(secondary.get(1));
            handler.enterChromosome("3");
            handler.secondary(secondary.get(2));
            handler.enterChromosome("4");
            handler.primary(primary.get(3));
        }};

        RegionZipper.zip(primary, secondary, handler);
    }

    private GenomeRegion createRegion(final String chromosome, final long start, final long end) {
        return ImmutableBEDGenomeRegion.builder().chromosome(chromosome).start(start).end(end).build();
    }
}
