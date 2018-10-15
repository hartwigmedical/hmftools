package com.hartwig.hmftools.common.variant.structural;

import static org.junit.Assert.assertEquals;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class StructuralVariantFactoryTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testSGL() {
        final String vcfEntry =
                "2\t192614842\tgridss14_291648b\tT\tTCTCTGTGACTTGAATGCAAACATCCCAAAGAAGTTTCTGAGAATGCTTCTGTCTAGATTTTCTCTGAAGACAATCCCGTTTCCAACGAAATCCTCAAGGCTAGGCAAATATCCTCTTGCAGATTCCAGAAAAAGAGTGTTTCAAAACTGCTCCTTCAAAACGGTGGTTCAATTCTCTTAGTTGAGTACACACATCTCAAATAAGTTTCTGAGAATGCTTCTGTCTAG.\t2076.59\tPASS\tAS=0;ASQ=0;ASRP=0;ASSR=0;BA=1;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=1002.16;BASRP=26;BASSR=14;BEALN=2:92296050|-|155M4I13M1I54M|6;BEID=asm14-297024;BEIDH=-1;BEIDL=0;BQ=2076.59;BSC=15;BSCQ=295.57;BUM=28;BUMQ=778.86;BVF=35;CAS=0;CASQ=0;CQ=2059.53;EVENT=gridss14_291648;IC=0;IQ=0;RAS=0;RASQ=0;REF=254;REFPAIR=88;RP=0;RPQ=0;SB=0.5862069;SC=1X;SR=0;SRQ=0;SVTYPE=BND;VF=0;BPI_AF=0.083743842364532;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0:0:0:0:0:0:0:0:0:0:27.37:0:0:1:27.37:1:0:0:0:0:0:43:14:0:0:0:0:0\t.:0:0:0:0:0:0:0:1002.16:26:14:2049.22:15:295.57:27:751.49:34:0:0:0:0:0:211:74:0:0:0:0:0";

        final StructuralVariant variant = StructuralVariantFactory.createSingleBreakend(codec.decode(vcfEntry));
        assertEquals(StructuralVariantType.SGL, variant.type());
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }
}
