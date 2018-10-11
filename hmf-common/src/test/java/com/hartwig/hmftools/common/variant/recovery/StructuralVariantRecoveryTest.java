package com.hartwig.hmftools.common.variant.recovery;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.region.GenomeRegion;
import com.hartwig.hmftools.common.region.GenomeRegionFactory;
import com.hartwig.hmftools.common.variant.structural.StructuralVariant;
import com.hartwig.hmftools.common.variant.structural.StructuralVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;

import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class StructuralVariantRecoveryTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testOrientation() {
        assertEquals(1, StructuralVariantRecovery.orientation("C[17:59493156["));
        assertEquals(-1, StructuralVariantRecovery.orientation("]17:59493156]C"));
    }

    @Test
    public void testMate() {
        assertEquals("17:59493156", StructuralVariantRecovery.mateLocation("C[17:59493156["));
        assertEquals("17:59493156", StructuralVariantRecovery.mateLocation("]17:59493156]C"));
    }

    @Test
    public void testClosestRegion() {
        final GenomeRegion start = GenomeRegionFactory.create("1", 1, 10000);
        final GenomeRegion middle1 = GenomeRegionFactory.create("1", 10001, 20000);
        final GenomeRegion middle2 = GenomeRegionFactory.create("1", 20001, 30000);
        final GenomeRegion end = GenomeRegionFactory.create("1", 30001, 40000);
        final List<GenomeRegion> regions = Lists.newArrayList(start, middle1, middle2, end);

        assertEquals(start, StructuralVariantRecovery.closest(1, regions));
        assertEquals(start, StructuralVariantRecovery.closest(5001, regions));
        assertEquals(middle1, StructuralVariantRecovery.closest(5002, regions));
        assertEquals(middle1, StructuralVariantRecovery.closest(15001, regions));
        assertEquals(middle2, StructuralVariantRecovery.closest(15002, regions));
        assertEquals(middle2, StructuralVariantRecovery.closest(25001, regions));
        assertEquals(end, StructuralVariantRecovery.closest(25002, regions));
        assertEquals(end, StructuralVariantRecovery.closest(40000, regions));
    }

    @Test
    public void testPete() {

        final String line1 =
                "16	33040203	gridss81_16816o	G	GGGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGC[16:33040204[	343.95	PON;BPI.Filter.SRSupportZero;qual	AS=0;ASQ=0;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=95.05;BSC=4;BSCQ=95.05;BUM=0;BUMQ=0;BVF=4;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816h;RAS=1;RASQ=343.95;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.6818182;SC=91M1X;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:136:107:0:0:0:0:0	.:0:5:11:0:0:0:0:0:0:0:95.05:4:95.05:0:0:4:0:0:0:343.95:343.95:941:672:0:0:0:0:16";
        final String line2 =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGCG	343.95	PON;BPI.Filter.SRSupportZero;qual	AS=1;ASQ=343.95;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=1240.5;BSC=57;BSCQ=1091.5;BUM=6;BUMQ=148.99;BVF=46;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816o;RAS=0;RASQ=0;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.5979381;SC=1X63M1D112M;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:60.92:2:35.02:1:25.9:3:0:0:0:0:0:136:107:0:0:0:0:0	.:343.95:5:11:0:0:0:0:0:0:0:1179.58:55:1056.49:5:123.09:43:0:0:0:343.95:0:941:672:0:0:0:0:16";

        final StructuralVariant variant = StructuralVariantFactory.create(codec.decode(line1), codec.decode(line2));
        System.out.println(variant);

    }

    @Test
    public void testPete2() {

        final String line1 =
                "11	437034	gridss60_16459o	A	AGCGGGGGGTGAGCGGGGGGTGAGCGGGGGGTGC[11:437035[	1313.53	PON	AS=0;ASQ=0;ASRP=0;ASSR=10;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm60-449095;BEIDH=0;BEIDL=0;BQ=637.58;BSC=6;BSCQ=130.5;BUM=19;BUMQ=507.08;BVF=21;CAS=0;CASQ=0;CQ=1276.16;EVENT=gridss60_16459;IC=18;IQ=843.44;PARID=gridss60_16459h;RAS=1;RASQ=432.72;REF=21;REFPAIR=40;RP=0;RPQ=0;SB=0.95;SC=109M1X;SR=2;SRQ=37.37;SVTYPE=BND;VF=18;BPI_AF=0.545454545454545,0.545454545454545;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:132.68:0:0:5:132.68:5:0:0:0:0:0:6:8:0:0:0:0:0	.:0:0:10:0:0:0:0:0:0:0:504.91:6:130.5:14:374.4:16:0:18:843.44:1313.53:432.72:15:32:0:0:2:37.37:18";
        final String line2 =
                "11	437035	gridss60_16459h	G	]11:437034]GCGGGGGGTGAGCGGGGGGTGAGCGGGGGGTGCG	1313.53	PON	AS=1;ASQ=432.72;ASRP=0;ASSR=10;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm60-449095;BEIDH=0;BEIDL=0;BQ=153.96;BSC=3;BSCQ=46.98;BUM=4;BUMQ=106.98;BVF=4;CAS=0;CASQ=0;CQ=1276.16;EVENT=gridss60_16459;IC=18;IQ=843.44;PARID=gridss60_16459o;RAS=0;RASQ=0;REF=21;REFPAIR=40;RP=0;RPQ=0;SB=0.9459459;SC=1X168M;SR=2;SRQ=37.37;SVTYPE=BND;VF=18;BPI_AF=0.545454545454545,0.545454545454545;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:53.07:0:0:2:53.07:2:0:0:0:0:0:6:8:0:0:0:0:0	.:432.72:0:10:0:0:0:0:0:0:0:100.89:3:46.98:2:53.91:2:0:18:843.44:1313.53:0:15:32:0:0:2:37.37:18";
        final StructuralVariant variant = StructuralVariantFactory.create(codec.decode(line1), codec.decode(line2));
        System.out.println(variant.type());

    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

}
