package com.hartwig.hmftools.store

import com.hartwig.hmftools.gripss.store.Breakend
import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.StructuralVariantContext
import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripss.store.HotspotStore
import junit.framework.Assert.assertTrue
import org.junit.Test

class HotspotStoreTest {
    private val defaultContigs = (1..22).map { it.toString() } + "X" + "Y" + "MT" + "M"
    private val contigComparator = ContigComparator(defaultContigs)

    @Test
    fun testPromiscuous() {
        val entry1 = Breakend.fromBed("7\t140481194\t140500161\t.\t.\t+")
        val store = HotspotStore(contigComparator, listOf(entry1), listOf())

        val line1 = "7\t129351321\tgridss545ff_157o\tA\tA]7:140488972]\t289.12\tLOW_QUAL\tAS=1;ASC=228M122D2X;ASQ=118.32;ASRP=9;ASSR=2;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm545-10992,asm550-3799;BEIDH=145,0;BEIDL=0,2;BMQ=60.00;BMQN=60.00;BMQX=60.00;BQ=14.87;BSC=1;BSCQ=14.87;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CIPOS=0,1;CIRPOS=-1,0;CQ=289.12;EVENT=gridss545ff_157;HOMLEN=1;HOMSEQ=C;IC=0;IHOMPOS=0,1;IQ=0.00;MATEID=gridss545ff_157h;MQ=58.13;MQN=50.00;MQX=60.00;RAS=1;RASQ=67.35;REF=92;REFPAIR=52;RP=5;RPQ=84.19;SB=0.5;SC=350M2X;SR=1;SRQ=19.26;SVTYPE=BND;VF=5\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:45:17:0:0.00:0:0.00:0\t.:118.32:9:2:0:0.00:0:0.00:0.00:0:0:14.87:1:14.87:0:0.00:0:0.00:0:0.00:289.12:67.35:47:35:5:84.19:1:19.26:5"
        val line2 = "7\t140488972\tgridss545ff_157h\tG\tG]7:129351321]\t289.12\tLOW_QUAL\tAS=1;ASC=325M1X;ASQ=67.35;ASRP=9;ASSR=2;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm545-10992,asm550-3799;BEIDH=0,2;BEIDL=145,0;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CIPOS=-1,0;CIRPOS=0,1;CQ=289.12;EVENT=gridss545ff_157;HOMLEN=1;HOMSEQ=G;IC=0;IHOMPOS=-1,0;IQ=0.00;MATEID=gridss545ff_157o;MQ=58.13;MQN=50.00;MQX=60.00;RAS=1;RASQ=118.32;REF=80;REFPAIR=40;RP=5;RPQ=84.19;SB=0.6666667;SC=325M1X;SR=1;SRQ=19.26;SVTYPE=BND;VF=5\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:31:13:0:0.00:0:0.00:0\t.:67.35:9:2:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:289.12:118.32:49:27:5:84.19:1:19.26:5"
        val context1 = VariantContextTestFactory.decode(line1)
        val context2 = VariantContextTestFactory.decode(line2)
        val sv1 = StructuralVariantContext(context1)
        val sv2 = StructuralVariantContext(context2)
        assertTrue(store.contains(sv1))
        assertTrue(store.contains(sv2))
    }

}