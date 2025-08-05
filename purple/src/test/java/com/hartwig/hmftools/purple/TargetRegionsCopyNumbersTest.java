package com.hartwig.hmftools.purple;

import static org.junit.Assert.assertEquals;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.cobalt.CobaltRatio;
import com.hartwig.hmftools.common.cobalt.CobaltRatioFile;
import com.hartwig.hmftools.common.genome.chromosome.Chromosome;
import com.hartwig.hmftools.common.purple.Gender;
import com.hartwig.hmftools.common.purple.PurpleCopyNumber;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;

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
        String panelFilePath =
                Paths.get("src", "test", "resources", "panel", "panel_1.bed").toAbsolutePath().toString();
        targetRegionsData.loadTargetRegionsBed(panelFilePath, ensemblDataCache);
        String cobaltDataFilePath = Resources.getResource("cobalt/sample_1.ratios.tsv").getPath();
        String purpleDataFilePath = Resources.getResource("purple/segments_1.tsv").getPath();

        Map<Chromosome, List<CobaltRatio>> cobaltData = CobaltRatioFile.readWithGender(cobaltDataFilePath, Gender.FEMALE, true);
        List<PurpleCopyNumber> purpleCopyNumbers = PurpleCopyNumberFile.read(purpleDataFilePath);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(targetRegionsData, cobaltData, purpleCopyNumbers, refGenomeVersion);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();
        assertEquals(4, copyNumberData.size());

        TargetRegionsCopyNumber cn0 = copyNumberData.get(0);
        assertEquals(26694001, cn0.mCobaltRatio().position());
        assertEquals(3, cn0.mOverlappingRegions().size());
        assertEquals(26694021, cn0.mOverlappingRegions().get(0).start());
        assertEquals("ARID1A_0", cn0.mOverlappingRegions().get(0).mTag);
        assertEquals("ARID1A_1", cn0.mOverlappingRegions().get(1).mTag);
        assertEquals("ARID1A_2", cn0.mOverlappingRegions().get(2).mTag);
        assertEquals(1, cn0.mPurpleCopyNumber().start());
        assertEquals(63461950, cn0.mPurpleCopyNumber().end());

        TargetRegionsCopyNumber cn2 = copyNumberData.get(2);
        assertEquals(55020001, cn2.mCobaltRatio().position());
        assertEquals(1, cn2.mOverlappingRegions().size());
        assertEquals(55020021, cn2.mOverlappingRegions().get(0).start());
        assertEquals("EGFR_1", cn2.mOverlappingRegions().get(0).mTag);
        assertEquals(38402474, cn2.mPurpleCopyNumber().start());
        assertEquals(59554330, cn2.mPurpleCopyNumber().end());

        TargetRegionsCopyNumber cn3 = copyNumberData.get(3);
        assertEquals(55021001, cn3.mCobaltRatio().position());
        assertEquals(1, cn3.mOverlappingRegions().size());
        assertEquals(55020021, cn3.mOverlappingRegions().get(0).start());
        assertEquals(55021520, cn3.mOverlappingRegions().get(0).end());
        assertEquals("EGFR_1", cn3.mOverlappingRegions().get(0).mTag);
    }

    @Test
    public void segmentsBreakCobaltWindows() throws IOException
    {
        // NOTE that if a Cobalt window gets broken into parts by Purple segments:
        // 1. We simply assign the Cobalt data to each new segment, with no attempt to re-calculate it per sub-window.
        // 2. We also assign each sub-window the panel regions for the entire Cobalt segment, with no
        //    attempt to assign these more accurately to the sub-windows.
        String panelFilePath =
                Paths.get("src", "test", "resources", "panel", "panel_2.bed").toAbsolutePath().toString();
        targetRegionsData.loadTargetRegionsBed(panelFilePath, ensemblDataCache);

        String cobaltDataFilePath = Resources.getResource("cobalt/sample_2.ratios.tsv").getPath();
        String purpleDataFilePath = Resources.getResource("purple/segments_2.tsv").getPath();

        Map<Chromosome, List<CobaltRatio>> cobaltData = CobaltRatioFile.readWithGender(cobaltDataFilePath, Gender.FEMALE, true);
        List<PurpleCopyNumber> purpleCopyNumbers = PurpleCopyNumberFile.read(purpleDataFilePath);
        TargetRegionsCopyNumbers trc = new TargetRegionsCopyNumbers(targetRegionsData, cobaltData, purpleCopyNumbers, refGenomeVersion);

        List<TargetRegionsCopyNumber> copyNumberData = trc.copyNumbersData();
        assertEquals(13, copyNumberData.size());

        //        12	11801001	-1	733.03	-1	0.9016	-1	-1	0.4732
        //        12	11801067	11801187	ETV6_UP_STREAM
        //        12	1	11804000	0.1036	7	0.5478	1.0000	TELOMERE	NONE	BAF_WEIGHTED	16	0.4816	1	1	0.0000	0.1036
        TargetRegionsCopyNumber cn0 = copyNumberData.get(0);
        assertEquals("12", cn0.mCobaltRatio().chromosome());
        assertEquals(11801001, cn0.mCobaltRatio().position());
        assertEquals(1, cn0.mOverlappingRegions().size());
        assertEquals(11801068, cn0.mOverlappingRegions().get(0).start());
        assertEquals("ETV6_UP_STREAM", cn0.mOverlappingRegions().get(0).mTag);
        assertEquals(1, cn0.mPurpleCopyNumber().start());
        assertEquals(11804000, cn0.mPurpleCopyNumber().end());

        //        12	11803001	-1	553.924	-1	0.8582	-1	-1	0.5317
        //        12	11803061	11803094	ETV6_CODING
        //        12	1	11804000	0.1036	7	0.5478	1.0000	TELOMERE	NONE	BAF_WEIGHTED	16	0.4816	1	1	0.0000	0.1036
        TargetRegionsCopyNumber cn1 = copyNumberData.get(1);
        assertEquals("12", cn1.mCobaltRatio().chromosome());
        assertEquals(11803001, cn1.mCobaltRatio().position());
        assertEquals(1, cn1.mOverlappingRegions().size());
        assertEquals(11803062, cn1.mOverlappingRegions().get(0).start());
        assertEquals(11803094, cn1.mOverlappingRegions().get(0).end());
        assertEquals("ETV6_CODING", cn1.mOverlappingRegions().get(0).mTag);
        assertEquals(1, cn1.mPurpleCopyNumber().start());
        assertEquals(11804000, cn1.mPurpleCopyNumber().end());

        //        12	11804001	-1	906.694	-1	1.0857	-1	-1	0.3804
        //        12	11804453	11804573	ETV6_INTRONIC_LONG
        //        12	11804001	36356693	2.2018	26	0.5217	0.5513	NONE	CENTROMERE	BAF_WEIGHTED	70	0.4028	11804001	11804001	0.9879	1.2139
        TargetRegionsCopyNumber cn2 = copyNumberData.get(2);
        assertEquals("12", cn2.mCobaltRatio().chromosome());
        assertEquals(11804001, cn2.mCobaltRatio().position());
        assertEquals(1, cn2.mOverlappingRegions().size());
        assertEquals(11804454, cn2.mOverlappingRegions().get(0).start());
        assertEquals(11804573, cn2.mOverlappingRegions().get(0).end());
        assertEquals("ETV6_INTRONIC_LONG", cn2.mOverlappingRegions().get(0).mTag);
        assertEquals(11804001, cn2.mPurpleCopyNumber().start());
        assertEquals(36356693, cn2.mPurpleCopyNumber().end());

        //        19	1610001	-1	1.652	-1	-1	-1	-1	0.5999
        //        19	1610705	1610755	TCF3_CODING
        //        19	1610001	1611050	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn3 = copyNumberData.get(3);
        assertEquals("19", cn3.mCobaltRatio().chromosome());
        assertEquals(1610001, cn3.mCobaltRatio().position());
        assertEquals(1, cn3.mOverlappingRegions().size());
        assertEquals(1610706, cn3.mOverlappingRegions().get(0).start());
        assertEquals(1610755, cn3.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn3.mOverlappingRegions().get(0).mTag);
        assertEquals(1610001, cn3.mPurpleCopyNumber().start());
        assertEquals(1611050, cn3.mPurpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1610001	1611050	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn4 = copyNumberData.get(4);
        assertEquals("19", cn4.mCobaltRatio().chromosome());
        assertEquals(1611001, cn4.mCobaltRatio().position());
        assertEquals(1, cn4.mOverlappingRegions().size());
        assertEquals(1611706, cn4.mOverlappingRegions().get(0).start());
        assertEquals(1611755, cn4.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn4.mOverlappingRegions().get(0).mTag);
        assertEquals(1610001, cn4.mPurpleCopyNumber().start());
        assertEquals(1611050, cn4.mPurpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611051	1611150	6.1096	1	0.5616	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn5 = copyNumberData.get(5);
        assertEquals("19", cn5.mCobaltRatio().chromosome());
        assertEquals(1611051, cn5.mCobaltRatio().position());
        assertEquals(1, cn5.mOverlappingRegions().size());
        assertEquals(1611706, cn5.mOverlappingRegions().get(0).start());
        assertEquals(1611755, cn5.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn5.mOverlappingRegions().get(0).mTag);
        assertEquals(1611051, cn5.mPurpleCopyNumber().start());
        assertEquals(1611150, cn5.mPurpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611151	1611350	8.1096	1	0.5636	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn6 = copyNumberData.get(6);
        assertEquals("19", cn6.mCobaltRatio().chromosome());
        assertEquals(1611151, cn6.mCobaltRatio().position());
        assertEquals(1, cn6.mOverlappingRegions().size());
        assertEquals(1611706, cn6.mOverlappingRegions().get(0).start());
        assertEquals(1611755, cn6.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn6.mOverlappingRegions().get(0).mTag);
        assertEquals(1611151, cn6.mPurpleCopyNumber().start());
        assertEquals(1611350, cn6.mPurpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611351	1611450	1.1096	1	0.5646	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn7 = copyNumberData.get(7);
        assertEquals(1611351, cn7.mCobaltRatio().position());
        assertEquals(1, cn7.mOverlappingRegions().size());
        assertEquals(1611706, cn7.mOverlappingRegions().get(0).start());
        assertEquals(1611755, cn7.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn7.mOverlappingRegions().get(0).mTag);
        assertEquals(1611351, cn7.mPurpleCopyNumber().start());
        assertEquals(1611450, cn7.mPurpleCopyNumber().end());

        //        19	1611001	-1	910.689	-1	0.8963	-1	-1	0.6122
        //        19	1611705	1611755	TCF3_CODING
        //        19	1611451	1613000	9.1096	1	0.5656	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn8 = copyNumberData.get(8);
        assertEquals(1611451, cn8.mCobaltRatio().position());
        assertEquals(1, cn8.mOverlappingRegions().size());
        assertEquals(1611706, cn8.mOverlappingRegions().get(0).start());
        assertEquals(1611755, cn8.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn8.mOverlappingRegions().get(0).mTag);
        assertEquals(1611451, cn8.mPurpleCopyNumber().start());
        assertEquals(1613000, cn8.mPurpleCopyNumber().end());

        //        19	1612001	-1	1760.197	-1	1.1675	-1	-1	0.6448
        //        19	1612705	1612755	TCF3_CODING
        //        19	1611451	1613000	9.1096	1	0.5656	0.7935	NONE	NONE	BAF_WEIGHTED	1	0.6448	1612001	1612001	1.2616	4.8479
        TargetRegionsCopyNumber cn9 = copyNumberData.get(9);
        assertEquals(1612001, cn9.mCobaltRatio().position());
        assertEquals(1, cn9.mOverlappingRegions().size());
        assertEquals(1612706, cn9.mOverlappingRegions().get(0).start());
        assertEquals(1612755, cn9.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn9.mOverlappingRegions().get(0).mTag);
        assertEquals(1611451, cn9.mPurpleCopyNumber().start());
        assertEquals(1613000, cn9.mPurpleCopyNumber().end());

        //        19	1613001	-1	800.883	-1	0.6991	-1	-1	0.5696
        //        19	1613705	1613755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn10 = copyNumberData.get(10);
        assertEquals(1613001, cn10.mCobaltRatio().position());
        assertEquals(1, cn10.mOverlappingRegions().size());
        assertEquals(1613706, cn10.mOverlappingRegions().get(0).start());
        assertEquals(1613755, cn10.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn10.mOverlappingRegions().get(0).mTag);
        assertEquals(1613001, cn10.mPurpleCopyNumber().start());
        assertEquals(26181781, cn10.mPurpleCopyNumber().end());

        //        19	1614001	-1	900.779	-1	0.7882	-1	-1	0.6444
        //        19	1614705	1614755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn11 = copyNumberData.get(11);
        assertEquals(1614001, cn11.mCobaltRatio().position());
        assertEquals(1, cn11.mOverlappingRegions().size());
        assertEquals(1614706, cn11.mOverlappingRegions().get(0).start());
        assertEquals(1614755, cn11.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn11.mOverlappingRegions().get(0).mTag);
        assertEquals(1613001, cn11.mPurpleCopyNumber().start());
        assertEquals(26181781, cn11.mPurpleCopyNumber().end());

        //        19	1615001	-1	1118.506	-1	0.9274	-1	-1	0.6566
        //        19	1615705	1615755	TCF3_CODING
        //        19	1613001	26181781	1.3304	49	0.5213	0.7656	NONE	CENTROMERE	BAF_WEIGHTED	117	0.5667	1613001	1613001	0.3118	1.0185
        TargetRegionsCopyNumber cn12 = copyNumberData.get(12);
        assertEquals(1615001, cn12.mCobaltRatio().position());
        assertEquals(1, cn12.mOverlappingRegions().size());
        assertEquals(1615706, cn12.mOverlappingRegions().get(0).start());
        assertEquals(1615755, cn12.mOverlappingRegions().get(0).end());
        assertEquals("TCF3_CODING", cn12.mOverlappingRegions().get(0).mTag);
        assertEquals(1613001, cn12.mPurpleCopyNumber().start());
        assertEquals(26181781, cn12.mPurpleCopyNumber().end());
    }
}
