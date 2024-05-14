package com.hartwig.hmftools.esvee.caller;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

import com.hartwig.hmftools.esvee.alignment.AlternativeAlignment;

import org.junit.Test;

public class RepeatMaskTest
{
    @Test
    public void testAlternativeAlignments()
    {
        // first BWA and VCF tag serialisation
        String bwaLocTagStr = "16,+24008715,44S28M46S,0;X,+133232624,44S27M47S,0;12,+54042138,37S35M46S,2;4,-84437081,46S25M47S,0";

        List<AlternativeAlignment> alignmentList = AlternativeAlignment.fromLocationTag(bwaLocTagStr);
        assertEquals(4, alignmentList.size());

        String vcfTag = AlternativeAlignment.toVcfTag(alignmentList);

        List<AlternativeAlignment> parsedAlignmentList = AlternativeAlignment.fromVcfTag(vcfTag);
        assertEquals(4, parsedAlignmentList.size());

        for(int i = 0; i < alignmentList.size(); ++i)
        {
            AlternativeAlignment alignment = alignmentList.get(i);
            AlternativeAlignment parsedAlignment = parsedAlignmentList.get(i);
            assertTrue(alignment.matches(parsedAlignment));
        }
    }

}
