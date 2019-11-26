package com.hartwig.hmftools.linx.misc;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.linx.neoepitope.InsertSeqData;
import com.hartwig.hmftools.linx.neoepitope.InsertSequenceAnalyser;

import org.junit.Test;

public class InsertSequenceTest
{
    @Test
    public void testMatchingLogic()
    {
        InsertSequenceAnalyser isAnalyser = new InsertSequenceAnalyser("", "", 5, 10, 2);

        final String sampleId = "SampleId1";
        final List<InsertSeqData> seqDataList = Lists.newArrayList();
        seqDataList.add(new InsertSeqData(sampleId, 1, "AAAACCCCGG"));
        seqDataList.add(new InsertSeqData(sampleId, 2, "TAAAACCCCG"));
        seqDataList.add(new InsertSeqData(sampleId, 3, "TTAAAACCCC"));
        seqDataList.add(new InsertSeqData(sampleId, 4, "GGAATTCCGG"));
        seqDataList.add(new InsertSeqData(sampleId, 5, "GAATTCCG"));
        seqDataList.add(new InsertSeqData(sampleId, 6, "CGAATTCCGGA"));

        isAnalyser.addSequenceData(seqDataList);

        isAnalyser.run();

        assertEquals(2, isAnalyser.getMatchSequences().size());
        assertEquals("AAAACCCC", isAnalyser.getMatchSequences().get(0));
        assertEquals("GAATTCCG", isAnalyser.getMatchSequences().get(1));

    }
}