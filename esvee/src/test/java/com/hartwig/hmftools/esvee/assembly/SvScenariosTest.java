package com.hartwig.hmftools.esvee.assembly;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.esvee.TestUtils.createJunctionReads;
import static com.hartwig.hmftools.esvee.TestUtils.loadRefGenomeBases;

import java.util.List;

import com.hartwig.hmftools.common.genome.region.Orientation;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.test.ReadIdGenerator;
import com.hartwig.hmftools.esvee.assembly.types.Junction;

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
        Junction posJunction = new Junction(CHR_1, posJunctionPos, FORWARD);
        Junction negJunction = new Junction(CHR_1, negJunctionPos, REVERSE);

        List<SAMRecord> junctionReads1 = createJunctionReads(
                mRefGenome, mReadIds.nextId(), 50, CHR_1, posJunctionPos, FORWARD, CHR_1, negJunctionPos, REVERSE, 50);


    }


}
