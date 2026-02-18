package com.hartwig.hmftools.orange.algo.isofox;

import static org.junit.Assert.assertEquals;

import java.io.IOException;

import com.google.common.io.Resources;

import org.junit.Ignore;
import org.junit.Test;

public class IsofoxDataLoaderTest
{
    private static final String ISOFOX_DATA_DIR = Resources.getResource("isofox").getPath();

    private static final String ISOFOX_SUMMARY_CSV = ISOFOX_DATA_DIR + "/sample.summary.csv";
    private static final String ISOFOX_GENE_DATA_CSV = ISOFOX_DATA_DIR + "/sample.gene_data.csv";
    private static final String ISOFOX_FUSION_CSV = ISOFOX_DATA_DIR + "/sample.pass_fusions.csv";
    private static final String ISOFOX_ALT_SPLICE_JUNCTION_CSV = ISOFOX_DATA_DIR + "/sample.alt_splice_junc.csv";

    @Ignore
    @Test
    public void canLoadIsofoxData() throws IOException
    {
        IsofoxData isofox = IsofoxDataLoader.load("sample", ISOFOX_DATA_DIR);

        assertEquals(2, isofox.geneExpressions().size());
        assertEquals(2, isofox.novelSpliceJunctions().size());
        assertEquals(2, isofox.fusions().size());
    }
}