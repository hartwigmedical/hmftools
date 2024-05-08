package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.common.test.MockRefGenome.generateRandomBases;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.buildFlags;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.TestUtils.BLACKLIST_LOCATIONS;
import static com.hartwig.hmftools.esvee.prep.TestUtils.HOTSPOT_CACHE;
import static com.hartwig.hmftools.esvee.prep.TestUtils.REGION_1;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.CANDIDATE_SUPPORT;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.JUNCTION;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.NO_SUPPORT;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertNotNull;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.hartwig.hmftools.common.region.BaseRegion;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadType;

import org.apache.commons.compress.utils.Lists;
import org.junit.Test;

public class JunctionsTest
{
    private static final String REF_BASES = generateRandomBases(500);

    private final ChrBaseRegion mPartitionRegion;
    private final JunctionTracker mJunctionTracker;

    public JunctionsTest()
    {
        mPartitionRegion = new ChrBaseRegion(CHR_1, 1, 5000);
        mJunctionTracker = new JunctionTracker(mPartitionRegion, new PrepConfig(1000), HOTSPOT_CACHE, BLACKLIST_LOCATIONS);
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
    public void testBlacklistRegions()
    {
        BLACKLIST_LOCATIONS.addRegion(CHR_1, new BaseRegion(500, 1500));

        JunctionTracker junctionTracker = new JunctionTracker(mPartitionRegion, new PrepConfig(1000), HOTSPOT_CACHE, BLACKLIST_LOCATIONS);

        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 100), "30S70M"));

        PrepRead read2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 820, REF_BASES.substring(20, 120), "100M"));

        read1.setReadType(JUNCTION);
        read2.setReadType(JUNCTION);
        junctionTracker.processRead(read1);
        junctionTracker.processRead(read2);

        PrepRead suppRead1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 800, REF_BASES.substring(0, 73), "3S70M"));

        suppRead1.setReadType(CANDIDATE_SUPPORT);
        junctionTracker.processRead(suppRead1);

        junctionTracker.assignJunctionFragmentsAndSupport();

        assertTrue(junctionTracker.junctions().isEmpty());
    }

    private void addDiscordantCandidate(
            final List<ReadGroup> discordantCandidates, final String readId, final String chr1, int pos1, final String chr2, int pos2)
    {
        PrepRead read = PrepRead.from(createSamRecord(readId, chr1, pos1, chr2, pos2, true, false, null));
        read.setReadType(CANDIDATE_SUPPORT);
        read.record().setMateNegativeStrandFlag(true);
        discordantCandidates.add(new ReadGroup(read));
    }

    @Test
    public void testDiscordantGroups()
    {
        // 5 fragments are required to support a discordant junction, unassigned to other junctions
        List<ReadGroup> discordantCandidates = Lists.newArrayList();

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 5200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 421, CHR_2, 5000); // unrelated
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 431, CHR_1, 5100);

        List<JunctionData> junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(0, junctions.size());

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 441, CHR_1, 5050);

        junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(2, junctions.size());
        assertTrue(junctions.get(0).JunctionGroups.isEmpty());
        assertEquals(540, junctions.get(0).Position);
        assertEquals(FORWARD, junctions.get(0).Orient);
        assertEquals(5050, junctions.get(1).Position);
        assertEquals(REVERSE, junctions.get(1).Orient);

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 540, CHR_1, 105000); // unrelated
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 550, CHR_1, 5300);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 560, CHR_1, 6000); // too far
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 570, CHR_1, 5000);

        // initially none because they are assigned to existing junctions
        junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(0, junctions.size());

        discordantCandidates.forEach(x -> x.clearJunctionPositions());
        junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(2, junctions.size());
        assertEquals(5, junctions.get(0).SupportingGroups.size());

        // a local DEL needs more support
        discordantCandidates.clear();
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(0, junctions.size());

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        junctions = DiscordantGroups.formDiscordantJunctions(REGION_1, discordantCandidates, DEFAULT_MAX_FRAGMENT_LENGTH);
        assertEquals(2, junctions.size());
    }

}
