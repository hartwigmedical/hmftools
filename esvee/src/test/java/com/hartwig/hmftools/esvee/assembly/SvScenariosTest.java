package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.TestUtils.loadRefGenomeBases;

import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

public class SvScenariosTest
{
    private final MockRefGenome mRefGenome;
    private final ReadIdGenerator mReadIds;

    public SvScenariosTest()
    {
        mRefGenome = new MockRefGenome(true);
        // loadRefGenomeBases(mRefGenome, "test_genome_01.csv");
        mReadIds = new ReadIdGenerator();
    }

    @Test
    public void testSimpleDeletion()
    {
        // minimum reads support
        Junction posJunction = new Junction(CHR_1, 60, POS_ORIENT);
        Junction negJunction = new Junction(CHR_1, 60, NEG_ORIENT);

        // Read read1 = createSamRecord("READ_01", 39, refBases.substring(9, 29), "20M");


    }


}
