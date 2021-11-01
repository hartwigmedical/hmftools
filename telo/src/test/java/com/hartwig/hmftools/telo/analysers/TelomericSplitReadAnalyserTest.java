package com.hartwig.hmftools.telo.analysers;

import static com.hartwig.hmftools.telo.TeloUtils.hasTelomericContent;

import com.hartwig.hmftools.telo.TeloUtils;
import com.hartwig.hmftools.telo.breakend.TelomericBreakEnd;
import com.hartwig.hmftools.telo.breakend.TelomericSplitReadAnalyser;

import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;
import static junit.framework.TestCase.assertEquals;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class TelomericSplitReadAnalyserTest
{
    @Test
    public void testFindTelomericBreakEnd()
    {
        TelomericSplitReadAnalyser analyser = new TelomericSplitReadAnalyser();
        analyser.setBoundaryZone(2);
        analyser.setMinSoftClipLength(10);

        SAMRecord record = new SAMRecord(null);
        record.setReadName("TestRead");
        record.setReferenceName("1");
        record.setAlignmentStart(1000);
        record.setCigarString("20S30M");
        record.setMateReferenceName("2");
        record.setMateAlignmentStart(200);
        record.setReadPairedFlag(true);
        record.setFirstOfPairFlag(true);

        // positive strand
        record.setReadNegativeStrandFlag(false);

        // first 18 is telomeric
        record.setReadString("TAACCCTAACCCTAACCC" + "AGCTAGCTAGCTACGTTGCAACCTGACGAA");

        TelomericBreakEnd tbe = analyser.findTelomericBreakEnd(record);

        assertEquals(TelomericBreakEnd.Type.LEFT_C_TELOMERIC, tbe.getType());
        assertEquals("1", tbe.getChromosome());
        assertEquals(999, tbe.getPosition());

        // try the negative strand
        record.setReadNegativeStrandFlag(true);
        record.setReadString(TeloUtils.reverseComplementSequence(record.getReadString()));

        tbe = analyser.findTelomericBreakEnd(record);

        assertEquals(TelomericBreakEnd.Type.LEFT_C_TELOMERIC, tbe.getType());
        assertEquals("1", tbe.getChromosome());
        assertEquals(999, tbe.getPosition());

        // right soft clip
        record.setCigarString("30M20S");
        // last 18 is telomeric
        record.setReadString("AGCTAGCTAGCTACGTTGCAACCTGACGAA" + "TAACCCTAACCCTAACCC");

        record.setReadNegativeStrandFlag(false);
        tbe = analyser.findTelomericBreakEnd(record);

        assertEquals(TelomericBreakEnd.Type.RIGHT_C_TELOMERIC, tbe.getType());
        assertEquals("1", tbe.getChromosome());
        assertEquals(1030, tbe.getPosition());

        record.setReadNegativeStrandFlag(true);
        record.setReadString(TeloUtils.reverseComplementSequence(record.getReadString()));
        tbe = analyser.findTelomericBreakEnd(record);

        assertEquals(TelomericBreakEnd.Type.RIGHT_C_TELOMERIC, tbe.getType());
        assertEquals("1", tbe.getChromosome());
        assertEquals(1030, tbe.getPosition());

        // test G telomeric
        record.setReadNegativeStrandFlag(false);
        record.setReadString("TTAGGGTTAGGGTTAGGG" + "AGCTAGCTAGCTACGTTGCAACCTGACGAA");
        record.setCigarString("20S30M");
        tbe = analyser.findTelomericBreakEnd(record);

        assertEquals(TelomericBreakEnd.Type.LEFT_G_TELOMERIC, tbe.getType());
        assertEquals("1", tbe.getChromosome());
        assertEquals(999, tbe.getPosition());
    }

    @Test
    public void testLikelyTelomeric()
    {
        String readBasesG = "AGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGGGTTAGTGTTAGGG";
        String readBasesC = "CTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCTAACCCAAACCCTAACCCAAACCCTAACCCTAACCCTAAC";

        assertTrue(TeloUtils.isLikelyGTelomeric(readBasesG));
        assertFalse(TeloUtils.isLikelyCTelomeric(readBasesG));

        assertFalse(TeloUtils.isLikelyGTelomeric(readBasesC));
        assertTrue(TeloUtils.isLikelyCTelomeric(readBasesC));

        //readBasesC = "CCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAACCAAACCCAAACCCTAACCCAAACCCACACCCCCACACCAACCCCCACCCCCACCACAACACCCCCCCCCCCCCCCCCCCCACC";

    }



}
