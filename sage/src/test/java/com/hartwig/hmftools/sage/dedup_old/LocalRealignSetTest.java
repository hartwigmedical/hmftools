package com.hartwig.hmftools.sage.dedup_old;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.variant.hotspot.ImmutableVariantHotspotImpl;
import com.hartwig.hmftools.common.variant.hotspot.VariantHotspot;
import com.hartwig.hmftools.sage.candidate.Candidate;
import com.hartwig.hmftools.sage.dedup_old.LocalRealignSet;
import com.hartwig.hmftools.sage.misc.CandidateSerializationTest;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class LocalRealignSetTest
{
    @Test
    public void testExample()
    {
        Candidate candidate1 = CandidateSerializationTest.decode(
                "10\t56022476\t.\tA\tAC\t698\tPASS\tLPS=13618;LRS=325;RC=TCACGGGGGC;RC_IDX=2;RC_LF=GTATGGGTTC;RC_NM=3;RC_RF=ATGGGAGCTG;TIER=LOW_CONFIDENCE;TNC=CAT\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t0/0:42,0:0.00:42:1620,0:43,0:0,0,0,0,42,42:0:0,0,0:0,0,0,0,1008,1008:43\t0/1:0,0:0.00:0:0,0:0,0:0,0,0,0,0,0:0:0,0,0:0,0,0,0,0,0:0\t0/1:65,32:0.327:98:3703,188:97,5:23,7,2,0,65,98:0:0,0,0:611,87,5,0,1518,2221:102");
        Candidate candidate2 = CandidateSerializationTest.decode(
                "10\t56022477\t.\tT\tG\t49\tmin_tumor_qual\tLPS=13618;RC=TCACGGGGGC;RC_IDX=4;RC_LF=GTATGGGTTC;RC_NM=3;RC_REPC=5;RC_REPS=G;RC_RF=ATGGGAGCTG;REP_C=4;REP_S=G;TIER=LOW_CONFIDENCE;TNC=ATG\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t1/1:0,0:0.00:41:0,0:0,0:0,0,0,0,0,41:0:0,0,0:0,0,0,0,0,1000:42\t0/1:0,0:0.00:0:0,0:0,0:0,0,0,0,0,0:0:0,0,0:0,0,0,0,0,0:0\t0/1:3,26:0.271:96:168,154:5,5:2,1,2,21,3,96:0:0,0,0:49,0,0,708,31,2350:99");
        Candidate candidate3 = CandidateSerializationTest.decode(
                "10\t56022477\t.\tT\tCG\t661\tPASS\tLPS=13618;LRS=325;MH=G;RC=CACGGGGGC;RC_IDX=2;RC_LF=TATGGGTTCT;RC_MH=G;RC_NM=2;RC_REPC=5;RC_REPS=G;RC_RF=ATGGGAGCTG;REP_C=4;REP_S=G;TIER=LOW_CONFIDENCE;TNC=ATG\tGT:AD:AF:DP:RABQ:RAD:RC_CNT:RC_IPC:RC_JIT:RC_QUAL:RDP\t1/1:0,0:0.00:42:0,0:0,0:0,0,0,0,0,42:0:0,0,0:0,0,0,0,0,1020:42\t0/1:0,0:0.00:0:0,0:0,0:0,0,0,0,0,0:0:0,0,0:0,0,0,0,0,0:0\t0/1:3,30:0.306:98:168,1020:5,26:21,7,0,2,3,98:0:0,0,0:565,96,0,49,19,2243:99");

        assertPhased(true, candidate1, candidate2);
        assertPhased(false, candidate1, candidate3);
        assertPhased(false, candidate2, candidate3);
    }

    private static void assertPhased(boolean expectedPhased, Candidate candidate1, Candidate candidate2)
    {
        int offset = LocalRealignSet.adjustedOffset(candidate1.variant(), candidate2.variant());
        assertEquals(expectedPhased, candidate1.readContext().phased(offset, candidate2.readContext()));
    }

    @Test
    public void testOffsetSNVs()
    {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "C", "T");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));
    }

    @Test
    public void testDelToTheRight()
    {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "CAT", "T");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));
    }

    @Test
    public void testInsToTheRight()
    {
        VariantHotspot left = create(100, "A", "T");
        VariantHotspot right = create(102, "C", "TAT");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));
    }

    @Test
    public void testDelToTheLeft()
    {
        VariantHotspot left = create(100, "AA", "A");
        VariantHotspot right = create(102, "G", "T");
        assertEquals(-1, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "AA", "A");
        right = create(103, "G", "T");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "AAC", "A");
        right = create(103, "G", "T");
        assertEquals(-1, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "AAC", "A");
        right = create(104, "G", "T");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));
    }

    @Test
    public void testInsToTheLeft()
    {
        VariantHotspot left = create(100, "A", "AA");
        VariantHotspot right = create(101, "G", "T");
        assertEquals(-2, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "A", "AA");
        right = create(102, "G", "T");
        assertEquals(-3, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "A", "AAC");
        right = create(101, "G", "T");
        assertEquals(-3, LocalRealignSet.adjustedOffset(left, right));

        left = create(100, "A", "AAC");
        right = create(102, "G", "T");
        assertEquals(-4, LocalRealignSet.adjustedOffset(left, right));
    }

    @Test
    public void testTwoInsertsAtSameLocation()
    {
        VariantHotspot left = create(100, "A", "AAA");
        VariantHotspot right = create(100, "A", "AA");
        assertEquals(0, LocalRealignSet.adjustedOffset(left, right));
        assertEquals(0, LocalRealignSet.adjustedOffset(right, left));
    }

    @Test
    public void testTwoDeletesAtSameLocation()
    {
        VariantHotspot left = create(100, "AAA", "A");
        VariantHotspot right = create(100, "AA", "A");
        assertEquals(0, LocalRealignSet.adjustedOffset(left, right));
        assertEquals(0, LocalRealignSet.adjustedOffset(right, left));
    }

    @NotNull
    static VariantHotspot create(int position, String ref, String alt)
    {
        return ImmutableVariantHotspotImpl.builder().chromosome("1").ref(ref).alt(alt).position(position).build();
    }
}
