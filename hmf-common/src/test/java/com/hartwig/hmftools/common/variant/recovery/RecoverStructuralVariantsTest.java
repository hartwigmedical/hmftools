package com.hartwig.hmftools.common.variant.recovery;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

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

public class RecoverStructuralVariantsTest {

    private static final String SAMPLE = "sample";

    private VCFCodec codec;

    @Before
    public void setup() {
        codec = createTestCodec();
    }

    @Test
    public void testMate() {
        assertEquals("17:59493156", RecoverStructuralVariants.mateLocation("C[17:59493156["));
        assertEquals("17:59493156", RecoverStructuralVariants.mateLocation("]17:59493156]C"));
    }

    @Test
    public void testClosestRegion() {
        final GenomeRegion start = GenomeRegionFactory.create("1", 1, 10000);
        final GenomeRegion middle1 = GenomeRegionFactory.create("1", 10001, 20000);
        final GenomeRegion middle2 = GenomeRegionFactory.create("1", 20001, 30000);
        final GenomeRegion end = GenomeRegionFactory.create("1", 30001, 40000);
        final List<GenomeRegion> regions = Lists.newArrayList(start, middle1, middle2, end);

        assertEquals(start, RecoverStructuralVariants.closest(1, regions));
        assertEquals(start, RecoverStructuralVariants.closest(5001, regions));
        assertEquals(middle1, RecoverStructuralVariants.closest(5002, regions));
        assertEquals(middle1, RecoverStructuralVariants.closest(15001, regions));
        assertEquals(middle2, RecoverStructuralVariants.closest(15002, regions));
        assertEquals(middle2, RecoverStructuralVariants.closest(25001, regions));
        assertEquals(end, RecoverStructuralVariants.closest(25002, regions));
        assertEquals(end, RecoverStructuralVariants.closest(40000, regions));
    }

    @Test
    public void testFilterFilter() {
        final String eligible =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGCG	343.95	BPI.Filter.SRSupportZero;qual	AS=1;ASQ=343.95;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=1240.5;BSC=57;BSCQ=1091.5;BUM=6;BUMQ=148.99;BVF=46;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816o;RAS=0;RASQ=0;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.5979381;SC=1X63M1D112M;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:60.92:2:35.02:1:25.9:3:0:0:0:0:0:136:107:0:0:0:0:0	.:343.95:5:11:0:0:0:0:0:0:0:1179.58:55:1056.49:5:123.09:43:0:0:0:343.95:0:941:672:0:0:0:0:16";
        final String passing =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGCG	343.95	.	AS=1;ASQ=343.95;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=1240.5;BSC=57;BSCQ=1091.5;BUM=6;BUMQ=148.99;BVF=46;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816o;RAS=0;RASQ=0;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.5979381;SC=1X63M1D112M;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:60.92:2:35.02:1:25.9:3:0:0:0:0:0:136:107:0:0:0:0:0	.:343.95:5:11:0:0:0:0:0:0:0:1179.58:55:1056.49:5:123.09:43:0:0:0:343.95:0:941:672:0:0:0:0:16";
        final String ponFiltered =
                "16	33040203	gridss81_16816o	G	GGGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGC[16:33040204[	343.95	PON	AS=0;ASQ=0;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=95.05;BSC=4;BSCQ=95.05;BUM=0;BUMQ=0;BVF=4;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816h;RAS=1;RASQ=343.95;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.6818182;SC=91M1X;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:136:107:0:0:0:0:0	.:0:5:11:0:0:0:0:0:0:0:95.05:4:95.05:0:0:4:0:0:0:343.95:343.95:941:672:0:0:0:0:16";
        final String afFiltered =
                "16	33040204	gridss81_16816h	G	]16:33040203]GGCGGCGGGGCAAGAAGCCGGGGCGGCGGGGTAAGAATCCGCG	343.95	af	AS=1;ASQ=343.95;ASRP=5;ASSR=11;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm81-367633;BEIDH=0;BEIDL=0;BQ=1240.5;BSC=57;BSCQ=1091.5;BUM=6;BUMQ=148.99;BVF=46;CAS=0;CASQ=0;CQ=411.71;EVENT=gridss81_16816;IC=0;IQ=0;PARID=gridss81_16816o;RAS=0;RASQ=0;REF=1077;REFPAIR=779;RP=0;RPQ=0;SB=0.5979381;SC=1X63M1D112M;SR=0;SRQ=0;SVTYPE=BND;VF=16;BPI_AF=0.0167189132706374,0.0167189132706374;LOCAL_LINKED_BY=;REMOTE_LINKED_BY=	GT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF	.:0:0:0:0:0:0:0:0:0:0:60.92:2:35.02:1:25.9:3:0:0:0:0:0:136:107:0:0:0:0:0	.:343.95:5:11:0:0:0:0:0:0:0:1179.58:55:1056.49:5:123.09:43:0:0:0:343.95:0:941:672:0:0:0:0:16";

        assertTrue(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(eligible)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(passing)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(ponFiltered)));
        assertFalse(RecoverStructuralVariants.isAppropriatelyFiltered(codec.decode(afFiltered)));
    }

    @NotNull
    private static VCFCodec createTestCodec() {
        VCFCodec codec = new VCFCodec();
        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
        return codec;
    }

}
