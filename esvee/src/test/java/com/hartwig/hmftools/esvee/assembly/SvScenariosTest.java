package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.esvee.TestUtils.createJunctionReads;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.TestUtils.loadRefGenomeBases;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.esvee.common.Junction;
import com.hartwig.hmftools.esvee.read.Read;

import org.junit.Test;

import htsjdk.samtools.SAMRecord;

public class SvScenariosTest
{
    private final MockRefGenome mRefGenome;
    private final ReadIdGenerator mReadIds;

    public SvScenariosTest()
    {
        mRefGenome = new MockRefGenome(true);
        loadRefGenomeBases(mRefGenome, "/test_genome_01.csv");
        mReadIds = new ReadIdGenerator();
    }

    @Test
    public void testSimpleDeletion()
    {
        // minimum reads support
        int posJunctionPos = 200;
        int negJunctionPos = 500;
        Junction posJunction = new Junction(CHR_1, posJunctionPos, POS_ORIENT);
        Junction negJunction = new Junction(CHR_1, negJunctionPos, NEG_ORIENT);

        List<SAMRecord> junctionReads1 = createJunctionReads(
                mRefGenome, mReadIds.nextId(), 50, CHR_1, posJunctionPos, POS_ORIENT, CHR_1, negJunctionPos, NEG_ORIENT, 50);


    }


}
