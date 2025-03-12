package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.ReferenceData.populateHlaTranscripts;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.lilac.GeneCache;

import org.junit.Test;

public class HlaGeneCacheTest
{
    @Test
    public void testClass1GeneCache()
    {
        Map<String,TranscriptData> hlaTranscriptMap = Maps.newHashMap();
        populateHlaTranscripts(hlaTranscriptMap, RefGenomeVersion.V37);

        assertEquals(3, hlaTranscriptMap.size());

        GeneCache geneCache = new GeneCache(hlaTranscriptMap);

        assertEquals(3, geneCache.GeneIds.size());

        // TO-DO: check all values match the constants for HLA class 1
    }
}
