package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.fusion.FusionCommon.POS_STRAND;
import static com.hartwig.hmftools.common.test.GeneTestUtils.TRANS_ID_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.createTransExons;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.common.gene.TranscriptData;
import com.hartwig.hmftools.common.test.GeneTestUtils;
import com.hartwig.hmftools.isofox.common.GeneCollection;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.novel.AltSpliceJunctionFinder;
import com.hartwig.hmftools.isofox.unmapped.UmrFinder;

import org.junit.Test;

public class UnmappedReadsTest
{
    @Test
    public void testValidReads()
    {
        String chromosome = CHR_1;
        String geneId = GENE_ID_1;

        GeneData geneData = new GeneData(geneId, geneId, chromosome, (byte) 1, 100, 1500, "");

        TranscriptData transData = createTransExons(
                geneId, TRANS_ID_1, POS_STRAND, new int[] {100, 300, 500}, 100, 150, 550, true, "");

        GeneReadData gene = new GeneReadData(geneData);

        gene.setTranscripts(Lists.newArrayList(transData));

        IsofoxConfig config = new IsofoxConfig();
        config.Functions.add(IsofoxFunction.UNMAPPED_READS);
        UmrFinder umrFinder = new UmrFinder(config, null);
        GeneCollection genes = new GeneCollection(0, Lists.newArrayList(gene));

        umrFinder.setGeneData(genes);



    }

    @Test
    public void testInvalidReads()
    {


    }

}
