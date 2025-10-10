package com.hartwig.hmftools.purple.targeted;

import static java.util.List.of;

import static com.hartwig.hmftools.common.cobalt.CobaltRatioFile.readWithGender;
import static com.hartwig.hmftools.common.purple.Gender.FEMALE;
import static com.hartwig.hmftools.common.purple.GermlineStatus.DIPLOID;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HET_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.HOM_DELETION;
import static com.hartwig.hmftools.common.purple.GermlineStatus.NOISE;
import static com.hartwig.hmftools.common.purple.GermlineStatus.UNKNOWN;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleSegment;
import com.hartwig.hmftools.common.region.TaggedRegion;

import org.junit.Before;
import org.junit.Test;

public class TargetRegionsCopyNumbersTest extends TargetRegionsTestBase
{
    @Before
    public void setup()
    {
        super.setup();
    }

    @Test
    public void multipleChromosomes() throws IOException
    {
        targetRegionsData.loadTargetRegionsBed(panelFilePath("panel_1.bed"), ensemblDataCache);
        Map<Chromosome, List<CobaltRatio>> cobaltData = readWithGender(cobaltFilePath("sample_1.ratios.tsv"), FEMALE, true);
        List<PurpleCopyNumber> purpleCopyNumbers = PurpleCopyNumberFile.read(purpleFilePath("segments_1.tsv"));
        PurpleSegment ps0 = ps("chr1", 1001, 26694000, DIPLOID);
        PurpleSegment ps1 = ps("chr1", 26694001, 250_000_000, HOM_DELETION);
        PurpleSegment ps2 = ps("chr7", 1, 55020000, HET_DELETION);
        PurpleSegment ps3 = ps("chr7", 55020001, 55035000, UNKNOWN);
        PurpleSegment ps4 = ps("chr7", 1, 55035000, NOISE);

        List<PurpleSegment> segments = of(ps0, ps1, ps2, ps3, ps4);
        TargetRegionsDataSource dataSource = new TargetRegionsDataSource(targetRegionsData, refGenomeVersion, segments);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(dataSource, cobaltData, purpleCopyNumbers);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();
        assertEquals(4, copyNumberData.size());

        TargetRegionDataBuilder targetRegionDataBuilder = new TargetRegionDataBuilder(
                refGenomeVersion, targetRegionsData.targetRegions(), segments, cobaltData, purpleCopyNumbers);

        targetRegionDataBuilder.buildTargetRegionData();
        List<TargetRegionsCopyNumber> results = targetRegionDataBuilder.targetRegionData();
        copyNumberData = results;

        TargetRegionsCopyNumber cn0 = copyNumberData.get(0);
        assertEquals("chr1", cn0.cobaltRatio().chromosome());
        assertEquals(26694001, cn0.cobaltRatio().position());
        assertEquals(3, cn0.overlappingRegions().size());
        assertEquals(26694021, cn0.overlappingRegions().get(0).start());
        assertEquals("ARID1A_0", cn0.overlappingRegions().get(0).mTag);
        assertEquals("ARID1A_1", cn0.overlappingRegions().get(1).mTag);
        assertEquals("ARID1A_2", cn0.overlappingRegions().get(2).mTag);
        assertEquals(1, cn0.purpleCopyNumber().start());
        assertEquals(63461950, cn0.purpleCopyNumber().end());
        assertEquals(HOM_DELETION, cn0.germlineStatus());

        TargetRegionsCopyNumber cn2 = copyNumberData.get(2);
        assertEquals("chr7", cn2.cobaltRatio().chromosome());
        assertEquals(55020001, cn2.cobaltRatio().position());
        assertEquals(1, cn2.overlappingRegions().size());
        assertEquals(55020021, cn2.overlappingRegions().get(0).start());
        assertEquals("EGFR_1", cn2.overlappingRegions().get(0).mTag);
        assertEquals(38402474, cn2.purpleCopyNumber().start());
        assertEquals(59554330, cn2.purpleCopyNumber().end());
        assertEquals(UNKNOWN, cn2.germlineStatus());

        TargetRegionsCopyNumber cn3 = copyNumberData.get(3);
        assertEquals(55021001, cn3.cobaltRatio().position());
        assertEquals(1, cn3.overlappingRegions().size());
        assertEquals(55020021, cn3.overlappingRegions().get(0).start());
        assertEquals(55021520, cn3.overlappingRegions().get(0).end());
        assertEquals("EGFR_1", cn3.overlappingRegions().get(0).mTag);
    }

    @Test
    public void segmentsBreakCobaltWindows() throws IOException
    {
        // NOTE that if a Cobalt window gets broken into parts by Purple segments:
        // 1. We simply assign the Cobalt data to each new segment, with no attempt to re-calculate it per sub-window.
        // 2. We also assign each sub-window the panel regions for the entire Cobalt segment, with no
        //    attempt to assign these more accurately to the sub-windows.
        targetRegionsData.loadTargetRegionsBed(panelFilePath("panel_2.bed"), ensemblDataCache);

        Map<Chromosome, List<CobaltRatio>> cobaltData = readWithGender(cobaltFilePath("sample_2.ratios.tsv"), FEMALE, true);
        List<PurpleCopyNumber> purpleCopyNumbers = PurpleCopyNumberFile.read(purpleFilePath("segments_2.tsv"));

        // these won't match
        PurpleSegment ps0 = ps("chr12", 1001, 26694000, DIPLOID);
        PurpleSegment ps1 = ps("chr19", 1001, 26694000, DIPLOID);
        List<PurpleSegment> segments = List.of(ps0, ps1);

        TargetRegionsDataSource dataSource = new TargetRegionsDataSource(targetRegionsData, refGenomeVersion, segments);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(dataSource, cobaltData, purpleCopyNumbers);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();

        TargetRegionDataBuilder targetRegionDataBuilder = new TargetRegionDataBuilder(
                refGenomeVersion, targetRegionsData.targetRegions(), segments, cobaltData, purpleCopyNumbers);

        targetRegionDataBuilder.buildTargetRegionData();
        List<TargetRegionsCopyNumber> results = targetRegionDataBuilder.targetRegionData();
        // copyNumberData = results;

        assertEquals(13, copyNumberData.size());

        //        12	11801001	-1	733.03	-1	0.9016	-1	-1	0.4732
        //        12	11801067	11801187	ETV6_UP_STREAM
        //        12	1	11804000	0.1036	7	0.5478	1.0000	TELOMERE	NONE	BAF_WEIGHTED	16	0.4816	1	1	0.0000	0.1036
        TargetRegionsCopyNumber cn0 = copyNumberData.get(0);
        assertEquals("chr12", cn0.cobaltRatio().chromosome());
        assertEquals(11801001, cn0.cobaltRatio().position());
        assertEquals(1, cn0.overlappingRegions().size());
        assertEquals(11801068, cn0.overlappingRegions().get(0).start());
        assertEquals("ETV6_UP_STREAM", cn0.overlappingRegions().get(0).mTag);
        assertEquals(1, cn0.purpleCopyNumber().start());
        assertEquals(11804000, cn0.purpleCopyNumber().end());

        //        12	11803001	-1	553.924	-1	0.8582	-1	-1	0.5317
        //        12	11803061	11803094	ETV6_CODING
        //        12	1	11804000	0.1036	7	0.5478	1.0000	TELOMERE	NONE	BAF_WEIGHTED	16	0.4816	1	1	0.0000	0.1036
        TargetRegionsCopyNumber cn1 = copyNumberData.get(1);
        assertEquals("chr12", cn1.cobaltRatio().chromosome());
        assertEquals(11803001, cn1.cobaltRatio().position());
        assertEquals(1, cn1.overlappingRegions().size());
        assertEquals(11803062, cn1.overlappingRegions().get(0).start());
        assertEquals(11803094, cn1.overlappingRegions().get(0).end());
        assertEquals("ETV6_CODING", cn1.overlappingRegions().get(0).mTag);
        assertEquals(1, cn1.purpleCopyNumber().start());
        assertEquals(11804000, cn1.purpleCopyNumber().end());

        //        12	11804001	-1	906.694	-1	1.0857	-1	-1	0.3804
        //        12	11804453	11804573	ETV6_INTRONIC_LONG
        //        12	11804001	36356693	2.2018	26	0.5217	0.5513	NONE	CENTROMERE	BAF_WEIGHTED	70	0.4028	11804001	11804001	0.9879	1.2139
        TargetRegionsCopyNumber cn2 = copyNumberData.get(2);
        assertEquals("chr12", cn2.cobaltRatio().chromosome());
        assertEquals(11804001, cn2.cobaltRatio().position());
        assertEquals(1, cn2.overlappingRegions().size());
        assertEquals(11804454, cn2.overlappingRegions().get(0).start());
        assertEquals(11804573, cn2.overlappingRegions().get(0).end());
        assertEquals("ETV6_INTRONIC_LONG", cn2.overlappingRegions().get(0).mTag);
        assertEquals(11804001, cn2.purpleCopyNumber().start());
        assertEquals(36356693, cn2.purpleCopyNumber().end());

        //        19	1610001	-1	1.652	-1	-1	-1	-1	0.5999
        //        19	1610705	1610755	TCF3_CODING
        //        19	1610001	1611050	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn3 = copyNumberData.get(3);
        assertEquals("chr19", cn3.cobaltRatio().chromosome());
        assertEquals(1610001, cn3.cobaltRatio().position());
        assertEquals(1, cn3.overlappingRegions().size());
        assertEquals(1610706, cn3.overlappingRegions().get(0).start());
        assertEquals(1610755, cn3.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn3.overlappingRegions().get(0).mTag);
        assertEquals(1610001, cn3.purpleCopyNumber().start());
        assertEquals(1611050, cn3.purpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1610001	1611050	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn4 = copyNumberData.get(4);
        assertEquals("chr19", cn4.cobaltRatio().chromosome());
        assertEquals(1611001, cn4.cobaltRatio().position());
        assertEquals(1, cn4.overlappingRegions().size());
        assertEquals(1611706, cn4.overlappingRegions().get(0).start());
        assertEquals(1611755, cn4.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn4.overlappingRegions().get(0).mTag);
        assertEquals(1610001, cn4.purpleCopyNumber().start());
        assertEquals(1611050, cn4.purpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611051	1611150	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn5 = copyNumberData.get(5);
        assertEquals("chr19", cn5.cobaltRatio().chromosome());
        assertEquals(1611051, cn5.cobaltRatio().position());
        assertEquals(1, cn5.overlappingRegions().size());
        assertEquals(1611706, cn5.overlappingRegions().get(0).start());
        assertEquals(1611755, cn5.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn5.overlappingRegions().get(0).mTag);
        assertEquals(1611051, cn5.purpleCopyNumber().start());
        assertEquals(1611150, cn5.purpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611151	1611350	8.1096	1	0.5636	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn6 = copyNumberData.get(6);
        assertEquals("chr19", cn6.cobaltRatio().chromosome());
        assertEquals(1611151, cn6.cobaltRatio().position());
        assertEquals(1, cn6.overlappingRegions().size());
        assertEquals(1611706, cn6.overlappingRegions().get(0).start());
        assertEquals(1611755, cn6.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn6.overlappingRegions().get(0).mTag);
        assertEquals(1611151, cn6.purpleCopyNumber().start());
        assertEquals(1611350, cn6.purpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611351	1611450	1.1096	1	0.5646	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn7 = copyNumberData.get(7);
        assertEquals(1611351, cn7.cobaltRatio().position());
        assertEquals(1, cn7.overlappingRegions().size());
        assertEquals(1611706, cn7.overlappingRegions().get(0).start());
        assertEquals(1611755, cn7.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn7.overlappingRegions().get(0).mTag);
        assertEquals(1611351, cn7.purpleCopyNumber().start());
        assertEquals(1611450, cn7.purpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611451	1613000	9.1096	1	0.5656	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn8 = copyNumberData.get(8);
        assertEquals(1611451, cn8.cobaltRatio().position());
        assertEquals(1, cn8.overlappingRegions().size());
        assertEquals(1611706, cn8.overlappingRegions().get(0).start());
        assertEquals(1611755, cn8.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn8.overlappingRegions().get(0).mTag);
        assertEquals(1611451, cn8.purpleCopyNumber().start());
        assertEquals(1613000, cn8.purpleCopyNumber().end());

        //        19	1612001	-1	1760.197	-1	1.1675	-1	-1	0.6448
        //        19	1612705	1612755	TCF3_CODING
        //        19	1611451	1613000	9.1096	1	0.5656	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn9 = copyNumberData.get(9);
        assertEquals(1612001, cn9.cobaltRatio().position());
        assertEquals(1, cn9.overlappingRegions().size());
        assertEquals(1612706, cn9.overlappingRegions().get(0).start());
        assertEquals(1612755, cn9.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn9.overlappingRegions().get(0).mTag);
        assertEquals(1611451, cn9.purpleCopyNumber().start());
        assertEquals(1613000, cn9.purpleCopyNumber().end());

        //        19	1613001	-1	800.883	-1	0.6991	-1	-1	0.5696
        //        19	1613705	1613755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn10 = copyNumberData.get(10);
        assertEquals(1613001, cn10.cobaltRatio().position());
        assertEquals(1, cn10.overlappingRegions().size());
        assertEquals(1613706, cn10.overlappingRegions().get(0).start());
        assertEquals(1613755, cn10.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn10.overlappingRegions().get(0).mTag);
        assertEquals(1613001, cn10.purpleCopyNumber().start());
        assertEquals(26181781, cn10.purpleCopyNumber().end());

        //        19	1614001	-1	900.779	-1	0.7882	-1	-1	0.6444
        //        19	1614705	1614755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn11 = copyNumberData.get(11);
        assertEquals(1614001, cn11.cobaltRatio().position());
        assertEquals(1, cn11.overlappingRegions().size());
        assertEquals(1614706, cn11.overlappingRegions().get(0).start());
        assertEquals(1614755, cn11.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn11.overlappingRegions().get(0).mTag);
        assertEquals(1613001, cn11.purpleCopyNumber().start());
        assertEquals(26181781, cn11.purpleCopyNumber().end());

        //        19	1615001	-1	1118.506	-1	0.9274	-1	-1	0.6566
        //        19	1615705	1615755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn12 = copyNumberData.get(12);
        assertEquals(1615001, cn12.cobaltRatio().position());
        assertEquals(1, cn12.overlappingRegions().size());
        assertEquals(1615706, cn12.overlappingRegions().get(0).start());
        assertEquals(1615755, cn12.overlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn12.overlappingRegions().get(0).mTag);
        assertEquals(1613001, cn12.purpleCopyNumber().start());
        assertEquals(26181781, cn12.purpleCopyNumber().end());
    }
}
