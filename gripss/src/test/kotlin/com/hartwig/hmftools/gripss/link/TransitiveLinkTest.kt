package com.hartwig.hmftools.gripss.link

import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripss.VariantContextTestFactory.createImpreciseVariant
import com.hartwig.hmftools.gripss.VariantContextTestFactory.createVariant
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import com.hartwig.hmftools.gripss.store.LinkStore
import com.hartwig.hmftools.gripss.store.VariantStore
import junit.framework.Assert.assertEquals
import junit.framework.Assert.assertTrue
import org.junit.Test

class TransitiveLinkTest {

    @Test
    fun testDUPvsINScase() {
        val v1 = VariantContextTestFactory.decode("4\t7873818\tgridss23_17823o\tA\tAAGAGGAAGAGTATTATTAGACCCAGTGATTTGC[4:7873819[\t4572.04\tPASS\tAS=1;ASQ=2268.61;ASRP=1;ASSR=176;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEALN=4:7873819|+|33M|60;BEID=asm23-312191,asm23-82937;BEIDH=0,0;BEIDL=0,0;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CQ=4572.04;EVENT=gridss23_17823;IC=11;IHOMPOS=0,0;IQ=586.76;LOCAL_LINKED_BY=dsb671;PARID=gridss23_17823h;RAS=1;RASQ=1716.68;REF=206;REFPAIR=123;REMOTE_LINKED_BY;RP=0;RPQ=0.00;SB=0.5205479;SC=114M1X;SR=0;SRQ=0.00;SVTYPE=BND;TAF=0.353;VF=104\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:33.56:0:2:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:33.56:0.00:19:25:0:0.00:0:0.00:2\t.:2235.05:1:174:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:11:586.76:4538.48:1716.68:187:98:0:0.00:0:0.00:102").toSv()
        val v2 = VariantContextTestFactory.decode("4\t7873819\tgridss23_17823h\tA\t]4:7873818]AGAGGAAGAGTATTATTAGACCCAGTGATTTGCA\t4504.37\tPASS\tAS=1;ASQ=1682.84;ASRP=1;ASSR=172;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEALN=4:7873819|+|33M|60;BEID=asm23-312191,asm23-82937;BEIDH=0,0;BEIDL=0,0;BQ=375.97;BSC=23;BSCQ=347.56;BUM=1;BUMQ=28.41;BVF=20;CAS=0;CASQ=0.00;CQ=4572.04;EVENT=gridss23_17823;IC=11;IHOMPOS=0,0;IQ=586.76;LOCAL_LINKED_BY;PARID=gridss23_17823o;RAS=1;RASQ=2234.77;REF=206;REFPAIR=123;REMOTE_LINKED_BY=dsb671;RP=0;RPQ=0.00;SB=0.5082645;SC=1X142M;SR=0;SRQ=0.00;SVTYPE=BND;TAF=0.348;VF=102\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:2:0:0.00:0:0.00:0.00:0:0:16.93:1:16.93:0:0.00:1:0.00:0:0.00:33.56:33.56:19:25:0:0.00:0:0.00:2\t.:1682.84:1:170:0:0.00:0:0.00:0.00:0:0:359.04:22:330.63:1:28.41:19:0.00:11:586.76:4470.81:2201.21:187:98:0:0.00:0:0.00:100").toSv()
        val v3 = VariantContextTestFactory.decode("4\t7873820\tgridss23_35386o\tG\t]4:7873852]G\t1184.62\tPASS\tAS=0;ASQ=0.00;ASRP=0;ASSR=0;BA=0;BANRP=0;BANRPQ=0.00;BANSR=63;BANSRQ=1184.62;BAQ=0.00;BASRP=0;BASSR=0;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CIPOS=-1,0;CIRPOS=-1,0;CQ=1184.62;EVENT=gridss23_35386;HOMLEN=1;HOMSEQ=A;IC=0;IHOMPOS=-1,0;IQ=0.00;LOCAL_LINKED_BY=dsb671;PARID=gridss23_35386h;RAS=0;RASQ=0.00;REF=206;REFPAIR=97;REMOTE_LINKED_BY;RP=0;RPQ=0.00;SB=0.52380955;SC=2X119M;SR=63;SRQ=1184.62;SVTYPE=BND;TAF=0.249;VF=62\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:19:24:0:0.00:0:0.00:0\t.:0.00:0:0:0:0.00:63:1184.62:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:1184.62:0.00:187:73:0:0.00:63:1184.62:62").toSv()
        val v4 = VariantContextTestFactory.decode( "4\t7873852\tgridss23_35386h\tA\tA[4:7873820[\t1184.62\tPASS\tAS=0;ASQ=0.00;ASRP=0;ASSR=0;BA=0;BANRP=0;BANRPQ=0.00;BANSR=63;BANSRQ=1184.62;BAQ=0.00;BASRP=0;BASSR=0;BQ=608.74;BSC=38;BSCQ=580.34;BUM=1;BUMQ=28.41;BVF=38;CAS=0;CASQ=0.00;CIPOS=-1,0;CIRPOS=-1,0;CQ=1184.62;EVENT=gridss23_35386;HOMLEN=1;HOMSEQ=A;IC=0;IHOMPOS=-1,0;IQ=0.00;LOCAL_LINKED_BY;PARID=gridss23_35386o;RAS=0;RASQ=0.00;REF=206;REFPAIR=97;REMOTE_LINKED_BY=dsb671;RP=0;RPQ=0.00;SB=0.56435645;SC=120M2X;SR=63;SRQ=1184.62;SVTYPE=BND;TAF=0.247;VF=62\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:33.56:2:33.56:0:0.00:2:0.00:0:0.00:0.00:0.00:17:25:0:0.00:0:0.00:0\t.:0.00:0:0:0:0.00:63:1184.62:0.00:0:0:575.18:36:546.78:1:28.41:36:0.00:0:0.00:1184.62:0.00:189:72:0:0.00:63:1184.62:62").toSv()

        val variantStore = VariantStore(listOf(v1,v2,v3,v4))
        val victim = TransitiveLink(LinkStore(listOf()), variantStore)
        val result = victim.transitiveLink(v1)
        assertTrue(result.isNotEmpty())
        assertEquals("gridss23_35386h<PAIR>gridss23_35386o", result[0].toString())

        // Check it is all symmetric
        assertTrue(victim.transitiveLink(v1).isNotEmpty())
        assertTrue(victim.transitiveLink(v2).isNotEmpty())
        assertTrue(victim.transitiveLink(v3).isNotEmpty())
        assertTrue(victim.transitiveLink(v4).isNotEmpty())

    }

    @Test
    fun transitiveLinksShouldNotBreakAssemblies() {

        val v1start = createImpreciseVariant("1", 1000, "id1s", "A", "A[1:2000[", 1000, "id1e").toSv()
        val v1end = createImpreciseVariant("1", 2000, "id1e", "A", "]1:1000]A", 1000, "id1s").toSv()

        val v2start = createVariant("1", 1000, "id2s", "A", "A[1:1400[", 1000, "id2e").toSv()
        val v2end = createVariant("1", 1400, "id2e", "A", "]1:1000]A", 1000, "id2s").toSv()
        val v3start = createVariant("1", 1600, "id3s", "A", "A[1:2000[", 1000, "id3e").toSv()
        val v3end = createVariant("1", 2000, "id3e", "A", "]1:1600]A", 1000, "id3s").toSv()
        val variantStore = VariantStore(listOf(v1start, v1end, v2start, v2end, v3start, v3end))

        assertTrue(TransitiveLink(LinkStore(listOf()), variantStore).transitiveLink(v1start).isNotEmpty())

        val v3link = Link("assemblyLink", "id3s", "randomOther", 0, 0)
        assertTrue(TransitiveLink(LinkStore(listOf(v3link)), variantStore).transitiveLink(v1start).isEmpty())
    }

}
