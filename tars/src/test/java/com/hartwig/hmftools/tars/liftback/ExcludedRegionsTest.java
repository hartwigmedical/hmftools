package com.hartwig.hmftools.tars.liftback;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.region.ChrBaseRegion;

import org.junit.Test;

import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;

// fragment-level exclusion contract: the PRIMARY placement decides the whole fragment. A primary inside an
// excluded zone drops the fragment; a supp/unmapped read inside it does not (only the primary is tested),
// so contaminating fragments are removed whole with no orphaned mates/supps left for dedup.
public class ExcludedRegionsTest
{
    private static final String CHR = "1";

    private static ExcludedRegions regions(final int start, final int end)
    {
        final Map<String, List<ChrBaseRegion>> map = new HashMap<>();
        map.put(CHR, new ArrayList<>(List.of(new ChrBaseRegion(CHR, start, end))));
        return new ExcludedRegions(map);
    }

    private static SAMRecord primary(final String contig, final int start, final int end)
    {
        final SAMRecord record = new SAMRecord(new SAMFileHeader());
        record.setReadName("read");
        record.setReferenceName(contig);
        record.setAlignmentStart(start);
        record.setCigarString((end - start + 1) + "M");
        return record;
    }

    private static SAMRecord supplementary(final String contig, final int start, final int end)
    {
        final SAMRecord record = primary(contig, start, end);
        record.setSupplementaryAlignmentFlag(true);
        return record;
    }

    @Test
    public void testFragmentExclusion()
    {
        // primary inside the region -> drop the fragment
        assertTrue(regions(1000, 2000).fragmentExcluded(List.of(primary(CHR, 1500, 1600))));

        // primary outside the region -> keep
        assertFalse(regions(1000, 2000).fragmentExcluded(List.of(primary(CHR, 2500, 2600))));

        // supp falls in the region but the primary does not -- the primary decides, so keep
        assertFalse(regions(1000, 2000).fragmentExcluded(
                List.of(primary(CHR, 5000, 5100), supplementary(CHR, 1500, 1600))));

        // primary on a different chromosome -> keep
        assertFalse(regions(1000, 2000).fragmentExcluded(List.of(primary("2", 1500, 1600))));

        // read end just touches the region start -- inclusive overlap must fire
        assertTrue(regions(2000, 3000).fragmentExcluded(List.of(primary(CHR, 1900, 2000))));
    }
}
