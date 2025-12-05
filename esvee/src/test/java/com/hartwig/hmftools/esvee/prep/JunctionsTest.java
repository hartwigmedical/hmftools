package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.bam.SamRecordUtils.CONSENSUS_READ_ATTRIBUTE;
import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.esvee.TestUtils.DEFAULT_MAP_QUAL;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.REF_BASES_400;
import static com.hartwig.hmftools.esvee.TestUtils.buildFlags;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.INVALID_JUNC_INDEX;
import static com.hartwig.hmftools.esvee.prep.JunctionUtils.findJunctionIndex;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEPTH_WINDOW_SIZE;
import static com.hartwig.hmftools.esvee.prep.TestUtils.HOTSPOT_CACHE;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.JUNCTION;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.NO_SUPPORT;

import static org.junit.Assert.assertNotEquals;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.bam.SupplementaryReadData;
import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.common.test.SamRecordTestUtils;
import com.hartwig.hmftools.esvee.common.ReadIdTrimmer;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadGroupStatus;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

import org.junit.Test;

public class JunctionsTest
{
    protected static final String REF_BASES = generateRandomBases(500);

    private final ChrBaseRegion mPartitionRegion;
    private final JunctionTracker mJunctionTracker;
    private final DepthTracker mDepthTracker;

    public JunctionsTest()
    {
        mPartitionRegion = new ChrBaseRegion(CHR_1, 1, 5000);

        mDepthTracker = new DepthTracker(new BaseRegion(mPartitionRegion.start(), mPartitionRegion.end()), DEPTH_WINDOW_SIZE);

        mJunctionTracker = new JunctionTracker(
                mPartitionRegion, new PrepConfig(1000), mDepthTracker, HOTSPOT_CACHE);
    }

    private void addRead(final PrepRead read, final ReadType readType)
    {
        read.setReadType(readType);
        mJunctionTracker.processRead(read);
    }

    @Test
    public void testBasicJunctions()
    {
        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 100), "30S70M"));

        PrepRead read2 = PrepRead.from(createSamRecord(
                read1.id(), CHR_1, 820, REF_BASES.substring(20, 120), "100M",
                buildFlags(false, true, false)));

        addRead(read1, JUNCTION);
        addRead(read2, NO_SUPPORT);

        PrepRead suppRead1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 73), "3S70M"));

        assertTrue(suppRead1.record().getMateNegativeStrandFlag());
        assertFalse(suppRead1.record().getReadNegativeStrandFlag());

        addRead(suppRead1, CANDIDATE_SUPPORT);

        PrepRead read3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 950, REF_BASES.substring(0, 100), "30S70M"));

        PrepRead read4 = PrepRead.from(createSamRecord(
                read3.id(), CHR_1, 980, REF_BASES.substring(20, 120), "100M",
                buildFlags(false, true, false)));

        addRead(read3, JUNCTION);
        addRead(read4, NO_SUPPORT);

        PrepRead suppRead2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 950, REF_BASES.substring(0, 73), "3S70M"));

        addRead(suppRead2, CANDIDATE_SUPPORT);

        PrepRead read5 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 950, REF_BASES.substring(20, 120), "100M"));

        PrepRead read6 = PrepRead.from(createSamRecord(
                read5.id(), CHR_1, 980, REF_BASES.substring(0, 100), "70M30S"));

        addRead(read5, NO_SUPPORT);
        addRead(read6, JUNCTION);

        PrepRead suppRead3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 990, REF_BASES.substring(0, 63), "60M3S"));

        addRead(suppRead3, CANDIDATE_SUPPORT);

        PrepRead read7 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1010, REF_BASES.substring(10, 90), "50M30S"));

        PrepRead read8 = PrepRead.from(createSamRecord(
                read7.id(), CHR_1, 1010, REF_BASES.substring(0, 50), "50M"));

        addRead(read7, JUNCTION);
        addRead(read8, NO_SUPPORT);

        PrepRead suppRead4 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 990, REF_BASES.substring(0, 73), "70M3S"));

        addRead(suppRead4, CANDIDATE_SUPPORT);

        mJunctionTracker.assignJunctionFragmentsAndSupport();

        assertEquals(4, mJunctionTracker.junctions().size());

        JunctionData junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 800).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(REVERSE, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(1, junctionData.exactSupportFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 950).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(REVERSE, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(1, junctionData.exactSupportFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 1049).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(FORWARD, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(1, junctionData.exactSupportFragmentCount());
        assertEquals(2, junctionData.supportingFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 1059).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(FORWARD, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(2, junctionData.exactSupportFragmentCount());
        assertEquals(0, junctionData.supportingFragmentCount());
    }

    @Test
    public void testCandidateOnlyJunctions()
    {
        mJunctionTracker.clear();

        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 100), "30S70M"));

        read1.record().setMappingQuality(0); // will cause the junction to be filtered

        addRead(read1, JUNCTION);

        PrepRead suppRead1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 810, CHR_3, 100, true, false, null));

        addRead(suppRead1, CANDIDATE_SUPPORT);

        PrepRead suppRead2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 73), "3S70M"));

        suppRead2.record().setMappingQuality(0);

        addRead(suppRead2, CANDIDATE_SUPPORT);

        PrepRead read3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 950, CHR_2, 100, true, false, null));

        addRead(read3, CANDIDATE_SUPPORT);

        mJunctionTracker.assignFragments();

        List<ReadGroup> junctionGroups = mJunctionTracker.formUniqueAssignedGroups();

        assertEquals(0, junctionGroups.size());

        List<ReadGroup> remoteCandidateGroups = mJunctionTracker.getRemoteCandidateReadGroups();
        assertEquals(4, remoteCandidateGroups.size());
    }

    @Test
    public void testInternalDeletes()
    {
        // initial delete is too short
        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 10, REF_BASES.substring(0, 80), "20M10D50M"));

        addRead(read1, JUNCTION);

        // then a simple one
        PrepRead read2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, REF_BASES.substring(0, 80), "20M40D20M"));

        addRead(read2, JUNCTION);

        // with supporting reads - first is too short as an indel
        PrepRead suppRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, REF_BASES.substring(0, 80), "20M20D20M"));

        addRead(suppRead, CANDIDATE_SUPPORT);

        suppRead = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 120, REF_BASES.substring(0, 80), "20M20D20M"));

        addRead(suppRead, CANDIDATE_SUPPORT);

        // and a more complicated one
        // 5S10M2D10M3I10M35D10M2S from base 210: 10-19 match, 20-21 del, 22-31 match, ignore insert, 32-41 match, 42-76 del, 77-86 match

        PrepRead read3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 210, REF_BASES.substring(0, 1), "5S10M2D10M3I10M35D10M2S"));

        addRead(read3, JUNCTION);

        mJunctionTracker.assignJunctionFragmentsAndSupport();

        assertEquals(4, mJunctionTracker.junctions().size());

        JunctionData junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 119).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(FORWARD, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(1, junctionData.exactSupportFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 160).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(REVERSE, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(1, junctionData.exactSupportFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 241).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(FORWARD, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(0, junctionData.exactSupportFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 277).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(REVERSE, junctionData.Orient);
        assertEquals(1, junctionData.junctionFragmentCount());
        assertEquals(0, junctionData.exactSupportFragmentCount());
    }

    @Test
    public void testInternalInserts()
    {
        // first is too short
        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 10, REF_BASES.substring(0, 70), "20M10I50M"));

        addRead(read1, JUNCTION);

        // then a simple one
        PrepRead read2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 100, REF_BASES.substring(0, 70), "20M40I50M"));

        addRead(read2, JUNCTION);

        // and a more complicated one

        PrepRead read3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 210, REF_BASES.substring(0, 100), "5S10M2D10M3I10M35I10M2S"));

        addRead(read3, JUNCTION);

        mJunctionTracker.assignJunctionFragmentsAndSupport();

        assertEquals(4, mJunctionTracker.junctions().size());
        assertEquals(119, mJunctionTracker.junctions().get(0).Position);
        assertEquals(120, mJunctionTracker.junctions().get(1).Position);

        assertEquals(241, mJunctionTracker.junctions().get(2).Position);
        assertEquals(242, mJunctionTracker.junctions().get(3).Position);
    }

    @Test
    public void testShortTemplatedInsertions()
    {
        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, REF_BASES_400.substring(0, 130), "32S66M32S"));

        PrepRead read2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 200, REF_BASES_400.substring(0, 130), "32S66M32S"));

        addRead(read1, JUNCTION);
        addRead(read2, JUNCTION);

        mJunctionTracker.assignFragments();

        assertEquals(2, mJunctionTracker.junctions().size());

        JunctionData junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 200).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(REVERSE, junctionData.Orient);
        assertEquals(2, junctionData.junctionFragmentCount());

        junctionData = mJunctionTracker.junctions().stream().filter(x -> x.Position == 265).findFirst().orElse(null);
        assertNotNull(junctionData);
        assertEquals(FORWARD, junctionData.Orient);
        assertEquals(2, junctionData.junctionFragmentCount());
    }

    @Test
    public void testPrimarySupplementaryDuplicates()
    {
        // primary and supplementary with matching coords and mates

        String readBases = REF_BASES.substring(0, 50);

        String lowerCigar = "20S30M";
        String upperCigar = "20M30S";
        String mateCigar = "50M";

        SupplementaryReadData suppDataUpper = new SupplementaryReadData(
                CHR_3, 1000, SupplementaryReadData.SUPP_POS_STRAND, upperCigar, DEFAULT_MAP_QUAL);

        PrepRead primary1 = PrepRead.from(SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1000, readBases, lowerCigar, CHR_2, 1000, false,
                false, suppDataUpper, true, mateCigar));
        primary1.setReadType(JUNCTION);

        SupplementaryReadData suppDataLower = new SupplementaryReadData(
                CHR_1, 1000, SupplementaryReadData.SUPP_POS_STRAND, lowerCigar, DEFAULT_MAP_QUAL);

        PrepRead supp1 = PrepRead.from(SamRecordTestUtils.createSamRecord(
                primary1.id(), CHR_3, 1000, readBases, upperCigar, CHR_2, 1000, false,
                true, suppDataLower, true, mateCigar));
        supp1.setReadType(CANDIDATE_SUPPORT);

        PrepRead supp2 = PrepRead.from(SamRecordTestUtils.createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 1000, readBases, lowerCigar, CHR_2, 1000, false,
                true, suppDataUpper, true, mateCigar));
        supp2.setReadType(JUNCTION);

        PrepRead primary2 = PrepRead.from(SamRecordTestUtils.createSamRecord(
                supp2.id(), CHR_3, 1000, readBases, upperCigar, CHR_2, 1000, false,
                false, suppDataLower, true, mateCigar));
        primary2.setReadType(CANDIDATE_SUPPORT);

        ReadIdTrimmer readIdTrimmer = new ReadIdTrimmer(true);

        Map<String,ReadGroup> readGroupsMap = Maps.newHashMap();

        ReadGroup readGroup1 = new ReadGroup(primary1, readIdTrimmer.trim(primary1.id()));
        readGroupsMap.put(readGroup1.id(), readGroup1);

        ReadGroup readGroup2 = new ReadGroup(supp2, readIdTrimmer.trim(supp2.id()));
        readGroupsMap.put(readGroup2.id(), readGroup2);

        JunctionUtils.markSupplementaryDuplicates(readGroupsMap, readIdTrimmer);
        assertNotEquals(ReadGroupStatus.DUPLICATE, readGroup1.groupStatus());
        assertEquals(ReadGroupStatus.DUPLICATE, readGroup2.groupStatus());

        // test that the upper reads find the same duplicate
        readGroup1 = new ReadGroup(supp1, readIdTrimmer.trim(supp1.id()));
        readGroupsMap.put(readGroup1.id(), readGroup1);

        readGroup2 = new ReadGroup(primary2, readIdTrimmer.trim(primary2.id()));
        readGroupsMap.put(readGroup2.id(), readGroup2);

        JunctionUtils.markSupplementaryDuplicates(readGroupsMap, readIdTrimmer);
        assertNotEquals(ReadGroupStatus.DUPLICATE, readGroup1.groupStatus());
        assertEquals(ReadGroupStatus.DUPLICATE, readGroup2.groupStatus());

        // repeat again, checking that consensus reads are favoured over non-consensus
        supp2.record().setAttribute(CONSENSUS_READ_ATTRIBUTE, "2:0");
        primary2.record().setAttribute(CONSENSUS_READ_ATTRIBUTE, "2:0");

        readGroup1 = new ReadGroup(primary1, readIdTrimmer.trim(primary1.id()));
        readGroupsMap.put(readGroup1.id(), readGroup1);

        readGroup2 = new ReadGroup(supp2, readIdTrimmer.trim(supp2.id()));
        readGroupsMap.put(readGroup2.id(), readGroup2);

        JunctionUtils.markSupplementaryDuplicates(readGroupsMap, readIdTrimmer);
        assertEquals(ReadGroupStatus.DUPLICATE, readGroup1.groupStatus());
        assertNotEquals(ReadGroupStatus.DUPLICATE, readGroup2.groupStatus());

        readGroup1 = new ReadGroup(supp1, readIdTrimmer.trim(supp1.id()));
        readGroupsMap.put(readGroup1.id(), readGroup1);

        readGroup2 = new ReadGroup(primary2, readIdTrimmer.trim(primary2.id()));
        readGroupsMap.put(readGroup2.id(), readGroup2);

        JunctionUtils.markSupplementaryDuplicates(readGroupsMap, readIdTrimmer);
        assertEquals(ReadGroupStatus.DUPLICATE, readGroup1.groupStatus());
        assertNotEquals(ReadGroupStatus.DUPLICATE, readGroup2.groupStatus());
    }

    @Test
    public void testJunctionDataCreateAndLookup()
    {
        // test junction-finding logic
        String readBases = REF_BASES.substring(0, 100);
        String leftCigar = "30S70M";

        PrepRead read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 100, readBases, leftCigar));
        mJunctionTracker.addJunctionData(read);

        // ensure junction matches prevent duplicate junctions
        mJunctionTracker.addJunctionData(read);

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 150, readBases, leftCigar));
        mJunctionTracker.addJunctionData(read);

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 200, readBases, leftCigar));
        mJunctionTracker.addJunctionData(read);

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 250, readBases, leftCigar));
        mJunctionTracker.addJunctionData(read);

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 300, readBases, leftCigar));
        mJunctionTracker.addJunctionData(read);

        // ensure junction matches prevent duplicate junctions
        mJunctionTracker.addJunctionData(read);

        List<JunctionData> junctions = mJunctionTracker.junctions();
        assertEquals(5, junctions.size());

        // add junctions at matching positions but opposite orientations
        String rightCigar = "50M30S";

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 101, readBases, rightCigar));
        mJunctionTracker.addJunctionData(read);

        assertEquals(6, junctions.size());

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 201, readBases, rightCigar));
        mJunctionTracker.addJunctionData(read);

        assertEquals(7, junctions.size());

        // ensure junction matches prevent duplicate junctions
        mJunctionTracker.addJunctionData(read);

        int juncIndex = findJunctionIndex(junctions, 50, 10);
        assertEquals(INVALID_JUNC_INDEX, juncIndex);

        juncIndex = findJunctionIndex(junctions, 100, 10);
        assertEquals(0, juncIndex);

        // using fast-search logic
        juncIndex = findJunctionIndex(junctions, 100, 0);
        assertEquals(0, juncIndex);

        juncIndex = findJunctionIndex(junctions, 149, 0);
        assertEquals(0, juncIndex);

        juncIndex = findJunctionIndex(junctions, 150, 0);
        assertEquals(1, juncIndex);

        juncIndex = findJunctionIndex(junctions, 151, 0);
        assertEquals(2, juncIndex);

        juncIndex = findJunctionIndex(junctions, 200, 0);
        assertEquals(3, juncIndex);

        juncIndex = findJunctionIndex(junctions, 225, 0);
        assertEquals(3, juncIndex);

        juncIndex = findJunctionIndex(junctions, 249, 0);
        assertEquals(3, juncIndex);

        juncIndex = findJunctionIndex(junctions, 250, 0);
        assertEquals(4, juncIndex);

        juncIndex = findJunctionIndex(junctions, 299, 0);
        assertEquals(5, juncIndex);

        juncIndex = findJunctionIndex(junctions, 300, 0);
        assertEquals(6, juncIndex);

        juncIndex = findJunctionIndex(junctions, 301, 0);
        assertEquals(6, juncIndex);

        // test adding discordant junctions - they will replace a split junction if have sufficiently more support

        // first matches an existing, does not replace since the split junction has sufficient support
        juncIndex = findJunctionIndex(junctions, 200, 0);

        JunctionData splitJunction = junctions.get(juncIndex);
        splitJunction.addJunctionReadGroup(new ReadGroup(read));
        splitJunction.addJunctionReadGroup(new ReadGroup(read));

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 151, readBases, "50M"));
        JunctionData discJunction = new JunctionData(200, REVERSE, read);
        discJunction.markDiscordantGroup();

        mJunctionTracker.addDiscordantJunction(discJunction);
        assertEquals(7, junctions.size());
        assertEquals(splitJunction, junctions.get(juncIndex));

        // check replacing
        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 251, readBases, "50M"));
        discJunction = new JunctionData(300, REVERSE, read);
        discJunction.markDiscordantGroup();

        juncIndex = findJunctionIndex(junctions, 300, 0);

        mJunctionTracker.addDiscordantJunction(discJunction);
        assertEquals(7, junctions.size());
        assertEquals(discJunction, junctions.get(juncIndex));

        // add a new discordant group
        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 221, readBases, "50M"));
        discJunction = new JunctionData(270, REVERSE, read);
        discJunction.markDiscordantGroup();

        mJunctionTracker.addDiscordantJunction(discJunction);
        assertEquals(8, junctions.size());
        juncIndex = findJunctionIndex(junctions, 270, 0);
        assertEquals(discJunction, junctions.get(juncIndex));
    }
}
