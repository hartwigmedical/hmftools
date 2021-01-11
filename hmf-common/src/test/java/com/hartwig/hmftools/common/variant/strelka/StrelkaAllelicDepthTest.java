package com.hartwig.hmftools.common.variant.strelka;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class StrelkaAllelicDepthTest {

    private static final String TUMOR = "tumor";
    private static final String NORMAL = "normal";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Lists.newArrayList(NORMAL, TUMOR));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

    @Test
    public void testSNP() {
        final String line = "1\t13\t.\tG\tA\t.\tPASS\tAC=0;AF=0.0;AN=0;NT=ref;QSS=13;QSS_NT=13;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1;"
                + "set=snvs\tGT:AU:CU:DP:FDP:GU:SDP:SUBDP:TU\t0/0:0,0:0,1:19:2:17,20:0:0:0,0\t0/1:7,8:0,1:53:7:39,53:0:0:0,0";

        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(line));
        assertTrue(victim.getGenotype(TUMOR).hasAD());
        assertEquals(39, victim.getGenotype(TUMOR).getAD()[0]);
        assertEquals(7, victim.getGenotype(TUMOR).getAD()[1]);

        assertTrue(victim.getGenotype(NORMAL).hasAD());
        assertEquals(17, victim.getGenotype(NORMAL).getAD()[0]);
        assertEquals(0, victim.getGenotype(NORMAL).getAD()[1]);
    }

    @Test
    public void testNonSNPNonIndel() {
        final String line =
                "1\t10352\t.\tT\t.\t.\tQSS_ref\tNT=ref;QSS=1;QSS_NT=1;SGT=TT->TT;SOMATIC;TQSS=1;TQSS_NT=1\tDP:FDP:SDP:SUBDP:AU:CU:GU:TU\t"
                        + "62:17:37:0:0,4:0,5:0,0:45,107\t173:53:93:0:6,26:0,8:0,2:114,341";

        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(line));
        assertFalse(victim.getGenotype(TUMOR).hasAD());
        assertFalse(victim.getGenotype(NORMAL).hasAD());
    }

    @Test
    public void testSNPWithNonEssentialAttributeMissing() {
        final String lineWithNonEssentialAttributesMissing =
                "1\t13\t.\tG\tA\t.\tPASS\tAC=0;AF=0.0;AN=0;NT=ref;QSS=13;QSS_NT=13;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs\t"
                        + "GT:AU:DP:FDP:GU:SDP:SUBDP:TU\t0/0:0,0:19:2:17,20:0:0:0,0\t0/1:7,8:53:7:39,53:0:0:0,0";

        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(lineWithNonEssentialAttributesMissing));
        assertTrue(victim.getGenotype(TUMOR).hasAD());
        assertEquals(39, victim.getGenotype(TUMOR).getAD()[0]);
        assertEquals(7, victim.getGenotype(TUMOR).getAD()[1]);

        assertTrue(victim.getGenotype(NORMAL).hasAD());
        assertEquals(17, victim.getGenotype(NORMAL).getAD()[0]);
        assertEquals(0, victim.getGenotype(NORMAL).getAD()[1]);
    }

    @Test
    public void testSNPWithEssentialAttributeMissing() {
        final String lineWithEssentialAttributesMissing =
                "1\t13\t.\tG\tA\t.\tPASS\tAC=0;AF=0.0;AN=0;NT=ref;QSS=13;QSS_NT=13;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs\t"
                        + "GT:CU:DP:FDP:GU:SDP:SUBDP:TU\t0/0:0,1:19:2:17,20:0:0:0,0\t0/1:0,1:53:7:39,53:0:0:0,0";

        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(lineWithEssentialAttributesMissing));
        assertFalse(victim.getGenotype(TUMOR).hasAD());
        assertFalse(victim.getGenotype(NORMAL).hasAD());
    }

    @Test
    public void testIndel() {
        final String line =
                "1\t50\t.\tC\tCAG\t.\tPASS\tAC=0;AF=0.0;AN=0;IC=1;IHP=5;NT=ref;QSI=15;QSI_NT=15;RC=0;RU=AG;SGT=ref->het;SOMATIC;TQSI=1;"
                        + "TQSI_NT=1;set=indels\tGT:DP:DP2:DP50:FDP50:SUBDP50:TAR:TIR:TOR\t0/0:30:30:28.39:2.89:0.0:29,40:0,0:1,2\t"
                        + "0/1:66:66:65.8:3.11:0.0:57,103:9,11:0,4";
        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(line));
        assertTrue(victim.getGenotype(TUMOR).hasAD());
        assertEquals(57, victim.getGenotype(TUMOR).getAD()[0]);
        assertEquals(9, victim.getGenotype(TUMOR).getAD()[1]);

        assertTrue(victim.getGenotype(NORMAL).hasAD());
        assertEquals(29, victim.getGenotype(NORMAL).getAD()[0]);
        assertEquals(0, victim.getGenotype(NORMAL).getAD()[1]);
    }

    @Test
    public void testKeepExistingAD() {
        final String line =
                "1\t13\t.\tG\tA\t.\tPASS\tAC=0;AF=0.0;AN=0;NT=ref;QSS=13;QSS_NT=13;SGT=GG->AG;SOMATIC;TQSS=1;TQSS_NT=1;set=snvs\t"
                        + "GT:AU:CU:DP:FDP:GU:SDP:SUBDP:TU:AD\t0/0:0,0:0,1:19:2:17,20:0:0:0,0:1,2\t0/1:7,8:0,1:53:7:39,53:0:0:0,0:3,4";

        final VariantContext victim = StrelkaAllelicDepth.enrich(codec.decode(line));
        assertTrue(victim.getGenotype(TUMOR).hasAD());
        assertEquals(3, victim.getGenotype(TUMOR).getAD()[0]);
        assertEquals(4, victim.getGenotype(TUMOR).getAD()[1]);

        assertTrue(victim.getGenotype(NORMAL).hasAD());
        assertEquals(1, victim.getGenotype(NORMAL).getAD()[0]);
        assertEquals(2, victim.getGenotype(NORMAL).getAD()[1]);
    }
}
