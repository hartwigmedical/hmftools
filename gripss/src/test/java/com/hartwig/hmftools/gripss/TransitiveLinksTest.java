package com.hartwig.hmftools.gripss;

import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.NEG_ORIENT;
import static com.hartwig.hmftools.common.utils.sv.SvCommonUtils.POS_ORIENT;
import static com.hartwig.hmftools.gripss.GripssTestUtils.CHR_1;
import static com.hartwig.hmftools.gripss.GripssTestUtils.buildLinkAttributes;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSgl;
import static com.hartwig.hmftools.gripss.GripssTestUtils.createSv;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_CIPOS;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_IMPRECISE;
import static com.hartwig.hmftools.gripss.common.VcfUtils.VT_QUAL;
import static com.hartwig.hmftools.gripss.links.DsbLinkFinder.findBreaks;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.gripss.common.Breakend;
import com.hartwig.hmftools.gripss.common.GenotypeIds;
import com.hartwig.hmftools.gripss.common.SvData;
import com.hartwig.hmftools.gripss.links.AlternatePath;
import com.hartwig.hmftools.gripss.links.AlternatePathFinder;
import com.hartwig.hmftools.gripss.links.AssemblyLinks;
import com.hartwig.hmftools.gripss.links.Link;
import com.hartwig.hmftools.gripss.links.LinkStore;
import com.hartwig.hmftools.gripss.links.TransitiveLinkFinder;

import org.junit.Test;

public class TransitiveLinksTest
{
    private final GripssTestApplication mGripss;

    public TransitiveLinksTest()
    {
        mGripss = new GripssTestApplication();
    }

    @Test
    public void testSortByQual()
    {
        Map<String, Object> tumorOverrides = Maps.newHashMap();
        tumorOverrides.put(VT_QUAL, 1);

        Map<String, Object> attributeOverrides = Maps.newHashMap();
        attributeOverrides.put(VT_CIPOS, new int[] {-10,10});
        attributeOverrides.put(VT_IMPRECISE, "true");

        SvData var1 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1000, 2000, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, attributeOverrides, null, tumorOverrides);

        tumorOverrides.put(VT_QUAL, 100);

        SvData var2 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1001, 2001, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, null, null, tumorOverrides);

        tumorOverrides.put(VT_QUAL, 1000);

        SvData var3 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1002, 2002, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, null, null, tumorOverrides);

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2, var3));

        TransitiveLinkFinder transLinkFinder = new TransitiveLinkFinder(dataCache, new LinkStore());
        List<Link> links = transLinkFinder.findTransitiveLinks(var1.breakendStart());
        assertEquals(1, links.size());
        assertEquals(var3.breakendStart(), links.get(0).breakendStart());
        assertEquals(var3.breakendEnd(), links.get(0).breakendEnd());

        //List<AlternatePath> alternatePaths = AlternatePathFinder.findPaths(dataCache, new LinkStore());
        //LinkStore transitiveLinkStore = AlternatePathFinder.createLinkStore(alternatePaths);
        // assertTrue(!transitiveLinkStore.getBreakendLinksMap().isEmpty());
    }

    @Test
    public void testDupVsInsert()
    {
        // an INS and DUP are linked as transitive candidates
        Map<String, Object> attributeOverrides = Maps.newHashMap();
        attributeOverrides.put(VT_CIPOS, new int[] { -1, 1 });

        SvData dup = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1002, 1024, NEG_ORIENT, POS_ORIENT, "",
                mGripss.mGenotypeIds, attributeOverrides, null, null);

        SvData ins = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1000, 1001, POS_ORIENT, NEG_ORIENT, "AGAGGAAGAGTATTATTAGACAG",
                mGripss.mGenotypeIds, attributeOverrides, null, null);

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(dup, ins));

        TransitiveLinkFinder transLinkFinder = new TransitiveLinkFinder(dataCache, new LinkStore());
        assertEquals(1, transLinkFinder.findTransitiveLinks(dup.breakendStart()).size());
        assertEquals(1, transLinkFinder.findTransitiveLinks(ins.breakendStart()).size());
        assertEquals(1, transLinkFinder.findTransitiveLinks(dup.breakendEnd()).size());
        assertEquals(1, transLinkFinder.findTransitiveLinks(ins.breakendEnd()).size());
    }

    @Test
    public void testTransitiveLinksShouldNotBreakAssemblies()
    {
        Map<String, Object> attributeOverrides = Maps.newHashMap();
        attributeOverrides.put(VT_IMPRECISE, "true");

        SvData var1 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1000, 2000, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, attributeOverrides, null);

        attributeOverrides = buildLinkAttributes("asm1", "1");

        SvData var2 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1000, 1400, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, null, attributeOverrides);

        SvData var3 = createSv(
                mGripss.IdGen.nextEventId(), CHR_1, CHR_1, 1600, 2000, POS_ORIENT, NEG_ORIENT, "",
                mGripss.mGenotypeIds, attributeOverrides, null);

        SvData var4 = createSgl(
                mGripss.IdGen.nextEventId(), CHR_1, 1500, POS_ORIENT, "",
                mGripss.mGenotypeIds, null, null, null);

        SvDataCache dataCache = new SvDataCache();

        GripssTestUtils.loadSvDataCache(dataCache, Lists.newArrayList(var1, var2, var3, var4));

        // first without any assembled links
        TransitiveLinkFinder transLinkFinder = new TransitiveLinkFinder(dataCache, new LinkStore());
        List<Link> links = transLinkFinder.findTransitiveLinks(var1.breakendStart());
        assertEquals(3, links.size());

        LinkStore assemblyLinks = AssemblyLinks.buildAssembledLinks(Lists.newArrayList(var1, var2, var3, var4));

        transLinkFinder = new TransitiveLinkFinder(dataCache, assemblyLinks);
        assertTrue(transLinkFinder.findTransitiveLinks(var1.breakendStart()).isEmpty());
    }
    /*

    @Test
    fun testMinTransitiveLinkDistance() {
        val v1start = createImpreciseVariant("1", 1000, "id1s", "A", "A[1:2000[", 1000, "id1e").toSv()
        val v1end = createImpreciseVariant("1", 2000, "id1e", "A", "]1:1000]A", 1000, "id1s").toSv()

        val v2start = createVariant("1", 1000, "id2s", "A", "A[1:1400[", 1000, "id2e").toSv()
        val v2end = createVariant("1", 1490, "id2e", "A", "]1:1000]A", 1000, "id2s").toSv()
        val v3start = createVariant("1", 1510, "id3s", "A", "A[1:2000[", 1000, "id3e").toSv()
        val v3end = createVariant("1", 2000, "id3e", "A", "]1:1600]A", 1000, "id3s").toSv()
        val variantStore = VariantStore(listOf(v1start, v1end, v2start, v2end, v3start, v3end))
        assertTrue(TransitiveLink(LinkStore(listOf()), variantStore).transitiveLink(v1start).isEmpty())
    }

    @Test
    fun testAllowALittleBuffer() {
        val v1 = VariantContextTestFactory.decode("1\t210565345\tgridss7_35037o\tT\t]1:210580840]T\t124.03\tLOW_QUAL;NO_ASSEMBLY\tAS=0;ASC=1X539N1X;ASQ=0.00;ASRP=0;ASSR=0;BA=0;BANRP=7;BANRPQ=124.03;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BQ=71.36;BSC=4;BSCQ=71.36;BUM=0;BUMQ=0.00;BVF=3;CAS=0;CASQ=0.00;CIPOS=-270,270;CIRPOS=-206,207;CQ=124.03;EVENT=gridss7_35037;IC=0;IMPRECISE;IQ=0.00;PARID=gridss7_35037h;RAS=0;RASQ=0.00;REF=162;REFPAIR=76;RP=7;RPQ=124.03;SB=0.5;SC=1X539N1X;SR=0;SRQ=0.00;SVTYPE=BND;VF=7\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:21:17:0:0.00:0:0.00:0\t.:0.00:0:0:7:124.03:0:0.00:0.00:0:0:71.36:4:71.36:0:0.00:3:0.00:0:0.00:124.03:0.00:141:59:7:124.03:0:0.00:7").toSv()
        val v2 = VariantContextTestFactory.decode("1\t210580840\tgridss7_35037h\tT\tT[1:210565345[\t124.03\tLOW_QUAL;NO_ASSEMBLY\tAS=0;ASC=1X412N1X;ASQ=0.00;ASRP=0;ASSR=0;BA=0;BANRP=7;BANRPQ=124.03;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BQ=0.00;BSC=0;BSCQ=0.00;BUM=0;BUMQ=0.00;BVF=0;CAS=0;CASQ=0.00;CIPOS=-206,207;CIRPOS=-270,270;CQ=124.03;EVENT=gridss7_35037;IC=0;IMPRECISE;IQ=0.00;PARID=gridss7_35037o;RAS=0;RASQ=0.00;REF=116;REFPAIR=58;RP=7;RPQ=124.03;SC=1X412N1X;SR=0;SRQ=0.00;SVTYPE=BND;VF=7\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:16:10:0:0.00:0:0.00:0\t.:0.00:0:0:7:124.03:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:124.03:0.00:100:48:7:124.03:0:0.00:7").toSv()

        val v3 = VariantContextTestFactory.decode("1\t210565616\tgridss7_35038o\tC\t]1:210567359]C\t1445.18\tPASS\tAS=1;ASC=1X125M;ASQ=384.43;ASRP=45;ASSR=22;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm7-336592,asm7-336596,asm7-6437,asm7-6438,asm7-6628;BEIDH=80,368,544,0,452;BEIDL=269,557,418,138,641;BQ=45.73;BSC=3;BSCQ=45.73;BUM=0;BUMQ=0.00;BVF=1;CAS=3;CASQ=407.52;CQ=1445.18;EVENT=gridss7_35038;IC=0;IHOMPOS=0,0;IQ=0.00;PARID=gridss7_35038h;RAS=1;RASQ=421.85;REF=188;REFPAIR=88;RP=3;RPQ=53.15;SB=0.47058824;SC=1X125M;SR=9;SRQ=178.23;SVTYPE=BND;VF=26\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:27:25:0:0.00:0:0.00:0\t.:384.43:45:22:0:0.00:0:0.00:0.00:0:0:45.73:3:45.73:0:0.00:1:407.52:0:0.00:1445.18:421.85:161:63:3:53.15:9:178.23:26").toSv()
        val v4 = VariantContextTestFactory.decode("1\t210567359\tgridss7_35038h\tG\tG[1:210565616[\t1445.18\tPASS\tAS=1;ASC=82M1D106M1X;ASQ=421.85;ASRP=45;ASSR=22;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm7-336592,asm7-336596,asm7-6437,asm7-6438,asm7-6628;BEIDH=269,557,418,138,641;BEIDL=80,368,544,0,452;BQ=31.00;BSC=2;BSCQ=31.00;BUM=0;BUMQ=0.00;BVF=0;CAS=3;CASQ=407.52;CQ=1445.18;EVENT=gridss7_35038;IC=0;IHOMPOS=0,0;IQ=0.00;PARID=gridss7_35038o;RAS=1;RASQ=384.43;REF=185;REFPAIR=95;RP=3;RPQ=53.15;SB=0.48809522;SC=82M1D106M1X;SR=9;SRQ=178.23;SVTYPE=BND;VF=26\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:32:26:0:0.00:0:0.00:0\t.:421.85:45:22:0:0.00:0:0.00:0.00:0:0:31.00:2:31.00:0:0.00:0:407.52:0:0.00:1445.18:384.43:153:69:3:53.15:9:178.23:26").toSv()

        val v5 = VariantContextTestFactory.decode("1\t210567170\tgridss7_35042o\tG\t]1:210580634]ATAGTAGAGTTGCTTGTACTTGG\t1204.33\tPASS\tAS=1;ASC=1X81M1D107M;ASQ=392.90;ASRP=33;ASSR=24;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm7-336592,asm7-336596,asm7-336597,asm7-6437,asm7-6628;BEIDH=0,0,0,755,0;BEIDL=80,368,490,544,452;BQ=48.46;BSC=3;BSCQ=48.46;BUM=0;BUMQ=0.00;BVF=0;CAS=3;CASQ=230.34;CQ=1204.33;EVENT=gridss7_35042;IC=0;IHOMPOS=0,0;IQ=0.00;PARID=gridss7_35042h;RAS=1;RASQ=431.97;REF=178;REFPAIR=102;RP=7;RPQ=124.03;SB=0.54320985;SC=1X81M1D107M;SR=1;SRQ=25.10;SVTYPE=BND;VF=25\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:0.00:0:0.00:0:0.00:0:0.00:0:0.00:0.00:0.00:29:21:0:0.00:0:0.00:0\t.:392.90:33:24:0:0.00:0:0.00:0.00:0:0:48.46:3:48.46:0:0.00:0:230.34:0:0.00:1204.33:431.97:149:81:7:124.03:1:25.10:25").toSv()
        val v6 = VariantContextTestFactory.decode("1\t210580634\tgridss7_35042h\tT\tTATAGTAGAGTTGCTTGTACTTG[1:210567170[\t1204.33\tPASS\tAS=1;ASC=84M2D383M1X;ASQ=431.97;ASRP=33;ASSR=24;BA=0;BANRP=0;BANRPQ=0.00;BANSR=0;BANSRQ=0.00;BAQ=0.00;BASRP=0;BASSR=0;BEID=asm7-336592,asm7-336596,asm7-336597,asm7-6437,asm7-6628;BEIDH=80,368,490,544,452;BEIDL=0,0,0,755,0;BQ=83.64;BSC=4;BSCQ=59.78;BUM=1;BUMQ=23.87;BVF=2;CAS=3;CASQ=230.34;CQ=1204.33;EVENT=gridss7_35042;IC=0;IHOMPOS=0,0;IQ=0.00;PARID=gridss7_35042o;RAS=1;RASQ=392.90;REF=126;REFPAIR=81;RP=7;RPQ=124.03;SB=0.5487805;SC=84M2D383M1X;SR=1;SRQ=25.10;SVTYPE=BND;VF=25\tGT:ASQ:ASRP:ASSR:BANRP:BANRPQ:BANSR:BANSRQ:BAQ:BASRP:BASSR:BQ:BSC:BSCQ:BUM:BUMQ:BVF:CASQ:IC:IQ:QUAL:RASQ:REF:REFPAIR:RP:RPQ:SR:SRQ:VF\t.:0.00:0:0:0:0.00:0:0.00:0.00:0:0:23.87:0:0.00:1:23.87:1:0.00:0:0.00:0.00:0.00:16:18:0:0.00:0:0.00:0\t.:431.97:33:24:0:0.00:0:0.00:0.00:0:0:59.78:4:59.78:0:0.00:1:230.34:0:0.00:1204.33:392.90:110:63:7:124.03:1:25.10:25").toSv()
        val variantStore = VariantStore(listOf(v1, v2, v3, v4, v5, v6).sortedWith(Comparator { x, y -> x.start.compareTo(y.start) }))

        val links = LinkStore(AssemblyLink(variantStore.selectAll()))
        assertTrue(TransitiveLink(links, variantStore).transitiveLink(v1).isNotEmpty())
    }

    @Test
    fun testBreathFirstSearch() {
        val v1start = createImpreciseVariant("1", 1000, "id1s", "A", "A[1:2000[", 1000, "id1e").toSv()
        val v1end = createImpreciseVariant("1", 2000, "id1e", "A", "]1:1000]A", 1000, "id1s").toSv()

        val lowQualMatchStart = createVariant("1", 1000, "lqs", "A", "A[1:2000[", 1, "lqe").toSv()
        val lowQualMatchEnd = createVariant("1", 2000, "lqe", "A", "]1:1000]A", 1, "lqs").toSv()

        val v2start = createVariant("1", 1000, "id2s", "A", "A[1:1400[", 1000, "id2e").toSv()
        val v2end = createVariant("1", 1400, "id2e", "A", "]1:1000]A", 1000, "id2s").toSv()
        val v3start = createVariant("1", 1600, "id3s", "A", "A[1:2000[", 1000, "id3e").toSv()
        val v3end = createVariant("1", 2000, "id3e", "A", "]1:1600]A", 1000, "id3s").toSv()
        val vSgl = createVariant("1", 1500, "idSgl", "A", "A.", 1000, listOf("PASS")).toSv()
        val variantStore = VariantStore(listOf(v1start, v1end, v2start, v2end, vSgl, v3start, v3end, lowQualMatchStart, lowQualMatchEnd))

        val result = TransitiveLink(LinkStore(listOf()), variantStore).transitiveLink(v1start)
        assertEquals(1, result.size)
    }

    @Test
    fun testHightestQualFirst() {
        val v1start = createImpreciseVariant("1", 1000, "id1s", "A", "A[1:2000[", 1000, "id1e").toSv()
        val v1end = createImpreciseVariant("1", 2000, "id1e", "A", "]1:1000]A", 1000, "id1s").toSv()

        val lowQualMatchStart = createVariant("1", 1000, "lqs", "A", "A[1:2000[", 1, "lqe").toSv()
        val lowQualMatchEnd = createVariant("1", 2000, "lqe", "A", "]1:1000]A", 1, "lqs").toSv()

        val highQualMatchStart = createVariant("1", 1000, "hqs", "A", "A[1:2000[", 10, "hqe").toSv()
        val highQualMatchEnd = createVariant("1", 2000, "hqe", "A", "]1:1000]A", 10, "hqs").toSv()
        val variantStore = VariantStore(listOf(v1start, v1end, lowQualMatchStart, lowQualMatchEnd, highQualMatchStart, highQualMatchEnd))

        val result = TransitiveLink(LinkStore(listOf()), variantStore).transitiveLink(v1start)
        assertEquals("hqs<PAIR>hqe", result[0].toString())
    }

    @Test
    fun testOneTransitiveLinkOnly() {
        val v1start = createImpreciseVariant("1", 1000, "id1s", "A", "A[1:5000[", 1000, "id1e").toSv()
        val v1end = createImpreciseVariant("1", 5000, "id1e", "A", "]1:1000]A", 1000, "id1s").toSv()

        val v2start = createVariant("1", 1000, "id2s", "A", "A[1:2000[", 1000, "id2e").toSv()
        val v2end = createVariant("1", 2000, "id2e", "A", "]1:1000]A", 1000, "id2s").toSv()
        val v3start = createVariant("1", 4000, "id3s", "A", "A[1:5000[", 1000, "id3e").toSv()
        val v3end = createVariant("1", 5000, "id3e", "A", "]1:4000]A", 1000, "id3s").toSv()

        val v4start = createVariant("1", 2500, "id4s", "A", "A[1:3500[", 1000, "id4e").toSv()
        val v4end = createVariant("1", 3500, "id4e", "A", "]1:2500]A", 1000, "id4s").toSv()

        val v5start = createVariant("1", 2500, "id5s", "A", "A[1:3500[", 1000, "id5e").toSv()
        val v5end = createVariant("1", 3500, "id5e", "A", "]1:2500]A", 1000, "id5s").toSv()

        val variantStoreV4 = variantStore(v1start, v1end, v2start, v2end, v4start, v4end, v3start, v3end)
        val variantStoreV5 = variantStore(v1start, v1end, v2start, v2end, v5start, v5end, v3start, v3end)
        val variantStoreAll = variantStore(v1start, v1end, v2start, v2end, v4start, v5start, v4end, v5end, v3start, v3end)

        assertTrue(TransitiveLink(LinkStore(listOf()), variantStoreV4).transitiveLink(v1start).isNotEmpty())
        assertTrue(TransitiveLink(LinkStore(listOf()), variantStoreV5).transitiveLink(v1start).isNotEmpty())
        assertTrue(TransitiveLink(LinkStore(listOf()), variantStoreAll).transitiveLink(v1start).isEmpty())
    }

     */



}
