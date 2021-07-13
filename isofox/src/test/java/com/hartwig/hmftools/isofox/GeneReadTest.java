package com.hartwig.hmftools.isofox;

import static com.hartwig.hmftools.common.gene.GeneTestUtils.createEnsemblGeneData;
import static com.hartwig.hmftools.isofox.TestUtils.CHR_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_ID_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_NAME_1;
import static com.hartwig.hmftools.isofox.TestUtils.GENE_START_1;
import static com.hartwig.hmftools.isofox.TestUtils.POS_STRAND;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.hartwig.hmftools.common.gene.GeneData;
import com.hartwig.hmftools.isofox.common.GeneReadData;
import com.hartwig.hmftools.isofox.common.RegionReadData;

import org.junit.Test;

public class GeneReadTest
{
    @Test
    public void testGeneRegions()
    {
        GeneData geneData = createEnsemblGeneData(GENE_ID_1, GENE_NAME_1, CHR_1, POS_STRAND, GENE_START_1, GENE_START_1 + 500);
        GeneReadData geneReadData = new GeneReadData(geneData);

        RegionReadData region1 = new RegionReadData(CHR_1, 1000, 1100);
        geneReadData.addExonRegion(region1);

        RegionReadData region2 = new RegionReadData(CHR_1, 1050, 1250);
        geneReadData.addExonRegion(region2);

        geneReadData.setHasUnsplicedRegions();
        assertFalse(geneReadData.hasUnsplicedRegions());

        RegionReadData region3 = new RegionReadData(CHR_1, 1200, 1400); // overlaps with 2nd not the 1st region
        geneReadData.addExonRegion(region3);

        geneReadData.setHasUnsplicedRegions();
        assertFalse(geneReadData.hasUnsplicedRegions());

        RegionReadData region4 = new RegionReadData(CHR_1, 1500, 1700);
        geneReadData.addExonRegion(region4);

        geneReadData.setHasUnsplicedRegions();
        assertTrue(geneReadData.hasUnsplicedRegions());

        RegionReadData region5 = new RegionReadData(CHR_1, 1600, 1800); // still leaves an unspliced region 1400-1500
        geneReadData.addExonRegion(region5);

        geneReadData.setHasUnsplicedRegions();
        assertTrue(geneReadData.hasUnsplicedRegions());
    }

}
