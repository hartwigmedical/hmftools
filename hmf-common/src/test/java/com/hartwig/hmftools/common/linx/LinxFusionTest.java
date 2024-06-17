package com.hartwig.hmftools.common.linx;

import static org.junit.Assert.assertNotNull;

import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.gene.TranscriptRegionType;

import org.junit.Test;

public class LinxFusionTest
{
    @Test
    public void canGenerateContextForAllTranscriptRegionTypesAndKnownTypes()
    {
        for(TranscriptRegionType regionType : TranscriptRegionType.values())
        {
            for(KnownFusionType knownType : KnownFusionType.values())
            {
                assertNotNull(LinxFusion.context(regionType, knownType, 0));
            }
        }
    }
}