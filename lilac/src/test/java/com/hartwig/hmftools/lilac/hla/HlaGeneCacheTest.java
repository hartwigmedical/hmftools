package com.hartwig.hmftools.lilac.hla;

import static com.hartwig.hmftools.lilac.ReferenceData.loadHlaTranscripts;

import static org.junit.Assert.assertEquals;

import java.util.Map;

import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.lilac.GeneCache;
import com.hartwig.hmftools.lilac.GeneSelector;

import org.junit.Test;

public class HlaGeneCacheTest
{
    @Test
    public void testClass1GeneCache()
    {
        Map<HlaGene, TranscriptData> hlaTranscriptMap = loadHlaTranscripts(RefGenomeVersion.V37, GeneSelector.MHC_CLASS_1);

        assertEquals(3, hlaTranscriptMap.size());

        GeneCache geneCache = new GeneCache(hlaTranscriptMap);

        assertEquals(3, geneCache.GeneNames.size());

        // TO-DO: check all values match the constants for HLA class 1
    }
}
