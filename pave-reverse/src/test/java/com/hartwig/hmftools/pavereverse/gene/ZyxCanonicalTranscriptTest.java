package com.hartwig.hmftools.pavereverse.gene;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.pavereverse.ReversePaveTestBase;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

public class ZyxCanonicalTranscriptTest extends ReversePaveTestBase
{
    private GeneTranscript geneTranscript;

    @Before
    public void setUp()
    {
        GeneData gene = reversePave.EnsemblCache.getGeneDataByName(zyx);
        TranscriptData transcriptData = reversePave.EnsemblCache.getTranscriptData(gene.GeneId, zyxCanonical);
        geneTranscript = new GeneTranscript(gene, transcriptData);
    }

    @Test
    public void totalTranslatedLength()
    {
        assertEquals(573*3, geneTranscript.totalTranslatedLength());
    }

    @Test
    public void absolutePositionOfTranslatedBase()
    {
        assertEquals(143_381_572, geneTranscript.absolutePositionOfTranslatedBase(1));
    }
}
