package com.hartwig.hmftools.cobalt.diploid;

import static junit.framework.TestCase.assertEquals;

import java.util.List;

import com.hartwig.hmftools.cobalt.Chromosome;
import com.hartwig.hmftools.cobalt.ChromosomePositionCodec;
import com.hartwig.hmftools.cobalt.CobaltColumns;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

import htsjdk.samtools.util.Locatable;

import tech.tablesaw.api.*;

public class DiploidRegionLoaderTest
{
    @Test
    public void testBuildRatios()
    {
        String chr1 = "1";
        String chr2 = "2";
        ChromosomePositionCodec chromosomePosCodec = new ChromosomePositionCodec();
        DiploidRegionLoader victim = new DiploidRegionLoader(chromosomePosCodec);
        victim.accept(locatable("1", 1001, 3000));
        victim.accept(locatable("1", 5001, 6000));
        victim.accept(locatable("2", 1001, 3000));

        Table result = victim.build();
        assertEquals(5, result.rowCount());
        StringColumn chrColumn = result.stringColumn(CobaltColumns.CHROMOSOME);
        assertReadRatio("1", 1001, result.where(chrColumn.isEqualTo(chr1)).row(0));
        assertReadRatio("1", 2001, result.where(chrColumn.isEqualTo(chr1)).row(1));
        assertReadRatio("1", 5001, result.where(chrColumn.isEqualTo(chr1)).row(2));
        assertReadRatio("2", 1001, result.where(chrColumn.isEqualTo(chr2)).row(0));
        assertReadRatio("2", 2001, result.where(chrColumn.isEqualTo(chr2)).row(1));
    }

    private void assertReadRatio(@NotNull String contig, long position, @NotNull Row victim)
    {
        assertEquals(contig, victim.getString(CobaltColumns.CHROMOSOME));
        assertEquals(position, victim.getInt(CobaltColumns.POSITION));
    }

    @NotNull
    private static Locatable locatable(String contig, int start, int end)
    {
        return new Locatable()
        {
            @Override
            public String getContig()
            {
                return contig;
            }

            @Override
            public int getStart()
            {
                return start;
            }

            @Override
            public int getEnd()
            {
                return end;
            }
        };
    }
}
