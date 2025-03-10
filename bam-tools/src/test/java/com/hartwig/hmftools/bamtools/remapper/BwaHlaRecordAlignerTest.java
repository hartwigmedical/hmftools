package com.hartwig.hmftools.bamtools.remapper;

import java.util.List;

import com.hartwig.hmftools.common.codon.Nucleotides;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;

import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import org.mockito.ArgumentCaptor;
import org.mockito.Mockito;

import htsjdk.samtools.SAMRecord;

public class BwaHlaRecordAlignerTest extends RemapperTestBase
{

    private PairAligner mAligner;
    private BwaHlaRecordAligner mBwaAligner;

    @Before
    public void before()
    {
        mAligner = Mockito.mock(PairAligner.class);
        mBwaAligner = new BwaHlaRecordAligner(mAligner, samFileHeader(), RefGenomeVersion.V38);
    }

    @Test
    public void alignPairTest()
    {
        List<BwaMemAlignment> alignmentsForRecord12 = List.of(
                bwa("145,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,5,31354375,-289")
        );
        List<BwaMemAlignment> alignmentsForRecord13 = List.of(
                bwa("97,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,5,31354513,289"),
                bwa("2161,1,32916486,32916522,67,103,0,0,36,36,67S36M48S,36,null,5,31354513,0"));
        ImmutablePair<List<BwaMemAlignment>, List<BwaMemAlignment>> alignments =
                ImmutablePair.of(alignmentsForRecord13, alignmentsForRecord12);
        Mockito.when(mAligner.alignSequences(Mockito.any(), Mockito.any())).thenReturn(alignments);

        SAMRecord record12 = records.get(12);
        SAMRecord record13 = records.get(13);

        List<SAMRecord> returned = mBwaAligner.alignPair(new RecordPair(record13, record12));
        Assert.assertEquals(3, returned.size());

        ArgumentCaptor<byte[]> captor1 = ArgumentCaptor.forClass(byte[].class);
        ArgumentCaptor<byte[]> captor2 = ArgumentCaptor.forClass(byte[].class);
        Mockito.verify(mAligner, Mockito.times(1)).alignSequences(captor1.capture(), captor2.capture());
        Assert.assertArrayEquals(Nucleotides.reverseComplementBases(record13.getReadBases()), captor1.getValue()); // Note that the pair is (record13, record12).
        Assert.assertArrayEquals(record12.getReadBases(), captor2.getValue());
    }
}
