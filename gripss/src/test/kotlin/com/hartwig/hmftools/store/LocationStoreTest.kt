package com.hartwig.hmftools.store

import com.hartwig.hmftools.bedpe.Breakend
import com.hartwig.hmftools.bedpe.Breakpoint
import com.hartwig.hmftools.gripss.ContigComparator
import com.hartwig.hmftools.gripss.VariantContextTestFactory
import com.hartwig.hmftools.gripss.VariantContextTestFactory.toSv
import com.hartwig.hmftools.gripss.store.LocationStore
import junit.framework.Assert.assertTrue
import org.junit.Assert
import org.junit.Ignore
import org.junit.Test

class LocationStoreTest {

    private val contigComparator = ContigComparator(null)

    @Test
    fun testContainSingle() {
        val entry1 = Breakend.fromBed("1\t665\t666\t.\t440\t+")
        val entry2 = Breakend.fromBed("1\t766\t767\t.\t440\t+")
        val store = LocationStore(contigComparator, listOf(entry1, entry2), listOf())
        Assert.assertTrue(store.contains(Breakend("1", 666, 666, 1)))
        Assert.assertFalse(store.contains(Breakend("1", 667, 667, 1)))
        Assert.assertTrue(store.contains(Breakend("1", 666, 666, 1)))
    }

    @Test
    fun testOutOfOrder() {
        val entry1 = Breakend.fromBed("1\t665\t666\t.\t440\t+")
        val entry2 = Breakend.fromBed("1\t766\t767\t.\t440\t+")
        val store = LocationStore(contigComparator, listOf(entry2, entry1), listOf())
        Assert.assertTrue(store.contains(Breakend("1", 666, 666, 1)))
        Assert.assertFalse(store.contains(Breakend("1", 667, 667, 1)))
        Assert.assertTrue(store.contains(Breakend("1", 666, 666, 1)))
    }

    @Test
    fun testGoBackFarEnough() {
        val entry1 = Breakpoint(Breakend("1", 224, 233, 1), Breakend("MT", 335, 335, 1))
        val entry2 = Breakpoint(Breakend("1", 229, 230, 1), Breakend("7", 335, 335, 1))
        val entry3 = Breakpoint(Breakend("1", 300, 300, 1), Breakend("MT", 335, 335, 1))
        val store = LocationStore(contigComparator, listOf(), listOf(entry2, entry1, entry3))

        Assert.assertTrue(store.contains(Breakpoint(Breakend("1", 230, 230, 1), entry1.endBreakend)))
        Assert.assertTrue(store.contains(Breakpoint(Breakend("1", 231, 231, 1), entry1.endBreakend)))
        Assert.assertTrue(store.contains(Breakpoint(Breakend("1", 230, 230, 1), entry1.endBreakend)))
    }

    @Test
    fun testBreakendsReversed() {
        val entry1 = Breakpoint(Breakend("1", 224, 233, 1), Breakend("MT", 335, 335, 1))
        val entry2 = Breakpoint(Breakend("1", 229, 230, 1), Breakend("7", 335, 335, 1))
        val entry3 = Breakpoint(Breakend("1", 300, 300, 1), Breakend("MT", 335, 335, 1))
        val store = LocationStore(contigComparator, listOf(), listOf(entry2, entry1, entry3))

        Assert.assertTrue(store.contains(Breakpoint(entry1.endBreakend, entry1.startBreakend)))
        Assert.assertTrue(store.contains(Breakpoint(entry2.startBreakend, entry1.endBreakend)))
    }


    @Test
    fun testSingleHotspot() {

        val contigComparator = ContigComparator(null)
        val entry1 = Breakpoint.fromBedpe("18\t23596577\t23681181\tX\t48104751\t48126879\tSS18-SSX1\t0\t-\t-\t18\tX", contigComparator)
        val entry2 = Breakpoint.fromBedpe("18\t23596577\t23681181\tX\t52725945\t52746239\tSS18-SSX2\t0\t-\t+\t18\tX", contigComparator)

        val variantString = "18\t23605969\tgridss86_482409b\tC\t.CTTTCCCATCATTTTGTGGGCCAGATGCTTCTGGCACTTCCTCCGAATCATTTCCTTCCTCTGCTGGCTTCTTGGGCATGATCTTTATAATGTGAAGATCACAGATAAATAGTATCAGTGACATATCTATAGTGCTTTTGAGCTTACAAAGGGTCTTCACATGCATTAGCTTATTCAATGTTCTCAACAACACTGGGAGAGTTACACAGGCCTAAATTAGGAGAAACCTGGGAGGGGAGGTTAGAAGGGAAAGGAATGGCCTAAGTGAATATGGTTTCCAGGGATAGAATGCTTATCTTCCCACTCTTTTAGGACTGACATTCTTGCAAACAGCAAAAATCTCCATGTAATTGAGAGTGTGATATACAGATGATTTGGAGAAGAGTAGCATTCTAAGAATTCACAAGGTCTACAAAGGGAAGAGGTTCTGTAAAATACAAGGGATCCCATATAAGCTTCTAGACAGCTGCTGGGAGAGTAAATGTAAAAACATAGGGAGGCGACAAAACACTGCTGGGAAAGATGGTGTGGGGAGATGAATACAGGGAAGGGAGAGGGAAAGAAATGGTTTGCCGAAATTAATCTAGGCAGCAAACAAAGCAGTACCAC\t5264.78\tPASS\tAS=0;ASC=1X;ASQ=0;ASRP=0;ASSR=0;BA=1;BANRP=0;BANRPQ=0;BANSR=0;BANSRQ=0;BAQ=2605.26;BASRP=79;BASSR=26;BEALN=X:52729547|+|608M|0,X:52786405|-|608M|,X:48209476|+|608M|;BEID=asm86-583807;BEIDH=-1;BEIDL=608;BPI_AF=0.36653;BQ=5264.78;BSC=26;BSCQ=461.82;BUM=81;BUMQ=2197.7;BVF=92;CAS=0;CASQ=0;CQ=4368.85;EVENT=gridss86_482409;IC=0;IQ=0;LOCAL_LINKED_BY=bebeins9;PURPLE_AF=0.417;PURPLE_CN=2;PURPLE_CN_CHANGE=0.874;PURPLE_PLOIDY=0.832;RAS=0;RASQ=0;REF=158;REFG=AAACACTGCACTGAAAATATA;REFPAIR=81;REMOTE_LINKED_BY=.;RP=0;RPQ=0;SB=0.42307693;SC=1X;SR=0;SRQ=0;SVTYPE=BND;VF=0;INSRMRT=MIRc;INSRMRC=SINE/MIR;INSRMRO=-;INSRMP=0.144736842105263\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:0:45:35:0:0:0:0:0\t.:0:0:0:0:0:0:0:2605.26:79:26:5264.78:26:461.82:81:2197.7:92:0:0:0:0:0:113:46:0:0:0:0:0\n"
        val variant = VariantContextTestFactory.decode(variantString).toSv()

        val store1 = LocationStore(contigComparator, listOf(), listOf(entry1))
        Assert.assertFalse(store1.contains(variant))

        val store2 = LocationStore(contigComparator, listOf(), listOf(entry2))
        Assert.assertTrue(store2.contains(variant))
    }

    @Ignore
    fun testStuff() {
        val breakEnds = Breakend.fromBedFile("/Users/jon/hmf/resources/gridss_pon_single_breakend.bed")
        val breakpoints = Breakpoint.fromBedpeFile("/Users/jon/hmf/resources/gridss_pon_breakpoint.bedpe", contigComparator)

        val ponStore = LocationStore(contigComparator, breakEnds, breakpoints)
        for (breakEnd in breakEnds) {
            assertTrue(ponStore.contains(breakEnd))
        }

        for (breakpoint in breakpoints) {
            assertTrue(ponStore.contains(breakpoint))
        }
    }

}