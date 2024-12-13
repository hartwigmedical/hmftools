package com.hartwig.hmftools.bamtools.remapper;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion.V38;

import static org.junit.Assert.assertEquals;

import org.broadinstitute.hellbender.utils.bwa.BwaMemAlignment;
import org.junit.Test;

public class PreferredAlignmentTest extends RemapperTestBase
{
    @Test
    public void xaTagAbsent() {
        BwaMemAlignment bwa = bwa("16,5,31354513,31354664,0,151,60,1,146,84,151M,111C39,null,-1,-1,0");
        BwaMemAlignment mate = bwa("16,5,31354375,31354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        PreferredAlignment preferred = new PreferredAlignment(bwa, mate, V38);
        assertEquals(preferred.getSamFlag(), bwa.getSamFlag());
        assertEquals(preferred.getRefId(), bwa.getRefId());
        assertEquals(preferred.getRefStart(), bwa.getRefStart() + 1);
        assertEquals(preferred.getMapQual(), bwa.getMapQual());
        assertEquals(preferred.getCigar(), bwa.getCigar());

        assertEquals(preferred.getMDTag(), bwa.getMDTag());
        assertEquals(preferred.getAlignerScore(), bwa.getAlignerScore());
        assertEquals(preferred.getSuboptimalScore(), bwa.getSuboptimalScore());
        assertEquals(preferred.getNMismatches(), bwa.getNMismatches());
    }

    @Test
    public void xaTagHasNoHlaAlternatives() {
        String noHlaXA = "chr6,+41354510,151M,6;chr6,-49942532,151M,10;chr6,-47889334,151M,10;chr6,-41360067,151M,10;";
        BwaMemAlignment bwa = bwa("16,5,41354513,41354664,0,151,60,1,146,84,151M,111C39", noHlaXA, "-1,-1,0");
        BwaMemAlignment mate = bwa("16,5,41354375,41354419,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        PreferredAlignment preferred = new PreferredAlignment(bwa, mate, V38);
        assertEquals(bwa.getRefStart() + 1, preferred.getRefStart());
    }

    @Test
    public void xaTagHasMultipleHlaAlternatives() {
        // chr6,+29944126,151M,0;chr6,-31270331,151M,6
        String multiHlaXA = "chr6,+32270330,151M,6;chr6,+29942932,151M,10;chr6,-27889334,151M,10;chr6,-32491142,151M,10;";
        BwaMemAlignment bwa = bwa("16,5,31354513,31354664,0,151,60,1,146,84,151M,111C39",  multiHlaXA, "-1,-1,0");
        BwaMemAlignment mate = bwa("16,5,29943400,29943551,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        PreferredAlignment preferred = new PreferredAlignment(bwa, mate, V38);
        assertEquals(29942932, preferred.getRefStart());
    }

    @Test
    public void xaTagHasAUniqueHlaAlternative() {
        String oneHlaXA = "chr6,+33268949,151M,6;chr6,+29942632,35S116M,10;chr6,-27945532,151M,10;chr6,-36360067,151M,10;";
        BwaMemAlignment bwa = bwa("16,5,31354513,31354664,0,151,60,1,146,84,151M,111C39", oneHlaXA, "-1,-1,0");
        BwaMemAlignment mate = bwa("16,5,2995000,2995151,0,44,60,1,39,20,44M107S,22C21,null,-1,-1,0");
        PreferredAlignment preferred = new PreferredAlignment(bwa, mate, V38);
        assertEquals(29942632, preferred.getRefStart());
        assertEquals(10, preferred.getMapQual());
        assertEquals("35S116M", preferred.getCigar());

        assertEquals(bwa.getMDTag(), preferred.getMDTag());
        assertEquals(bwa.getAlignerScore(), preferred.getAlignerScore());
        assertEquals(bwa.getSuboptimalScore(), preferred.getSuboptimalScore());
        assertEquals(bwa.getNMismatches(), preferred.getNMismatches());
    }

    @Test
    public void baseAlignmentIsBetterThanAnyAlternative()
    {
        String xa = "chr6,-31356387,151M,8;chr6,+29888372,151M,9;";
        BwaMemAlignment bwa = bwa("81,5,31271292,31271443,0,151,18,6,121,111,151M,20G10G12G1G34A53C15", xa,"5,31270953,-490");
        BwaMemAlignment mate = bwa("161,5,31270953,31271104,0,151,60,7,119,83,151M,28C2T51T36T18T4G4C1,null,5,31271292,490");
        PreferredAlignment preferred = new PreferredAlignment(bwa, mate, V38);
        assertEquals(preferred.getSamFlag(), bwa.getSamFlag());
        assertEquals(preferred.getRefId(), bwa.getRefId());
        assertEquals(preferred.getRefStart(), bwa.getRefStart() + 1);
        assertEquals(preferred.getMapQual(), bwa.getMapQual());
        assertEquals(preferred.getCigar(), bwa.getCigar());
    }
}
