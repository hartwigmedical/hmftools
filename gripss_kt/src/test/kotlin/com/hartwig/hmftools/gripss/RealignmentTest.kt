package com.hartwig.hmftools.gripss

import com.hartwig.hmftools.gripss.VariantContextTestFactory.decode
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import junit.framework.Assert.assertEquals
import org.junit.Test

class RealignmentTest {
    private val defaultContigs = (1..22).map { it.toString() } + "X" + "Y" + "MT" + "M"
    private val contigComparator = ContigComparator(defaultContigs)

    @Test
    fun testSameOrientationIsSymmetric() {
        val start = decode("9\t96821078\tgridss54_33293o\tC\t[10:62178441[C\t2243.88\tPASS\tAS=1;ASC=2X501M;ASQ=754.14;ASRP=62;ASSR=31;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm54-279690,asm58-155360;BEIDH=0,503;BEIDL=555,0;BQ=169.48;BSC=5;BSCQ=82.82;BUM=3;BUMQ=86.66;BVF=2;CAS=0;CASQ=0;CIPOS=0,1;CIRPOS=-1,0;CQ=2227.89;EVENT=gridss54_33293;HOMLEN=1;HOMSEQ=C;IC=0;IHOMPOS=0,1;IQ=0;LOCAL_LINKED_BY;PARID=gridss54_33293h;RAS=1;RASQ=818.43;REF=150;REFPAIR=98;REMOTE_LINKED_BY;RP=32;RPQ=511.78;SB=0.613636;SC=2X501M;SR=8;SRQ=159.54;SVTYPE=BND;TAF=0.214;VF=48\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0:0:0:0:0:0:0:0:0:0:28.94:0:0:1:28.94:1:0:0:0:0:0:45:27:0:0:0:0:0\t.:754.14:62:31:0:0:0:0:0:0:0:140.53:5:82.82:2:57.72:1:0:0:0:2243.88:818.43:105:71:32:511.78:8:159.54:48").toSv()
        val end = decode("10\t62178441\tgridss54_33293h\tT\t[9:96821078[T\t2243.88\tPASS\tAS=1;ASC=2X553M;ASQ=818.43;ASRP=62;ASSR=31;BA=0;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=0;BASRP=0;BASSR=0;BEID=asm54-279690,asm58-155360;BEIDH=555,0;BEIDL=0,503;BQ=179.09;BSC=10;BSCQ=179.09;BUM=0;BUMQ=0;BVF=0;CAS=0;CASQ=0;CIPOS=-1,0;CIRPOS=0,1;CQ=2227.89;EVENT=gridss54_33293;HOMLEN=1;HOMSEQ=G;IC=0;IHOMPOS=-1,0;IQ=0;LOCAL_LINKED_BY;PARID=gridss54_33293o;RAS=1;RASQ=754.14;REF=80;REFPAIR=53;REMOTE_LINKED_BY;RP=32;RPQ=511.78;SB=0.489796;SC=2X553M;SR=8;SRQ=159.54;SVTYPE=BND;TAF=0.397;VF=48\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:35:25:0:0:0:0:0\t.:818.43:62:31:0:0:0:0:0:0:0:179.09:10:179.09:0:0:0:0:0:0:2243.88:754.14:45:28:32:511.78:8:159.54:48").toSv()
        val startRealigned = start.centreAlignConfidenceIntervals(contigComparator)
        val endRealigned = end.centreAlignConfidenceIntervals(contigComparator)

        assertEquals(startRealigned.first, endRealigned.second)
        assertEquals(startRealigned.second, endRealigned.first)
    }
}