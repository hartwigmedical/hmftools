package com.hartwig.hmftools.esvee.prep;

import static com.hartwig.hmftools.common.genome.region.Orientation.FORWARD;
import static com.hartwig.hmftools.common.genome.region.Orientation.REVERSE;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_1;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_2;
import static com.hartwig.hmftools.common.test.GeneTestUtils.CHR_3;
import static com.hartwig.hmftools.esvee.TestUtils.READ_ID_GENERATOR;
import static com.hartwig.hmftools.esvee.TestUtils.createSamRecord;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DEFAULT_MAX_FRAGMENT_LENGTH;
import static com.hartwig.hmftools.esvee.prep.PrepConstants.DISCORDANT_GROUP_MAX_LOCAL_LENGTH;
import static com.hartwig.hmftools.esvee.prep.TestUtils.REGION_1;
import static com.hartwig.hmftools.esvee.prep.types.ReadType.CANDIDATE_SUPPORT;

import static junit.framework.TestCase.assertEquals;
import static junit.framework.TestCase.assertFalse;
import static junit.framework.TestCase.assertTrue;

import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.region.ChrBaseRegion;
import com.hartwig.hmftools.esvee.prep.types.JunctionData;
import com.hartwig.hmftools.esvee.prep.types.PrepRead;
import com.hartwig.hmftools.esvee.prep.types.ReadGroup;

import org.junit.Test;

public class DiscordantGroupsTest
{
    private void addDiscordantCandidate(
            final List<ReadGroup> discordantCandidates, final String readId, final String chr1, int pos1, final String chr2, int pos2)
    {
        PrepRead read = PrepRead.from(createSamRecord(readId, chr1, pos1, chr2, pos2, true, false, null));
        read.setReadType(CANDIDATE_SUPPORT);
        read.record().setMateNegativeStrandFlag(true);
        discordantCandidates.add(new ReadGroup(read));
    }

    @Test
    public void testDiscordantReadCriteria()
    {
        List<KnownHotspot> knownHotspots = Lists.newArrayList();

        knownHotspots.add(new KnownHotspot(
                new ChrBaseRegion(CHR_1, 100, 1000), FORWARD, new ChrBaseRegion(CHR_2, 100, 1000), REVERSE,
                ""));

        DiscordantGroups discordantGroups = new DiscordantGroups(
                REGION_1, DEFAULT_MAX_FRAGMENT_LENGTH, knownHotspots, false);

        PrepRead read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 500, CHR_2, 500,
                true, false, null));

        ReadGroup readGroup = new ReadGroup(read);
        assertTrue(discordantGroups.isDiscordantGroup(readGroup));
        assertTrue(discordantGroups.isRelevantDiscordantGroup(readGroup));

        // outside the known pair ranges
        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 500, CHR_2, 10000,
                true, false, null));

        readGroup = new ReadGroup(read);
        assertTrue(discordantGroups.isDiscordantGroup(readGroup));
        assertFalse(discordantGroups.isRelevantDiscordantGroup(readGroup));

        // otherwise local within the required distance
        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 500, CHR_1, 10000,
                true, false, null));
        read.record().setInferredInsertSize(9500);

        readGroup = new ReadGroup(read);
        assertTrue(discordantGroups.isDiscordantGroup(readGroup));
        assertTrue(discordantGroups.isRelevantDiscordantGroup(readGroup));

        read = PrepRead.from(createSamRecord(READ_ID_GENERATOR.nextId(), CHR_1, 500, CHR_1, DISCORDANT_GROUP_MAX_LOCAL_LENGTH * 2,
                true, false, null));
        read.record().setInferredInsertSize(read.record().getMateAlignmentStart() - read.start());

        readGroup = new ReadGroup(read);
        assertTrue(discordantGroups.isDiscordantGroup(readGroup));
        assertFalse(discordantGroups.isRelevantDiscordantGroup(readGroup));
    }

    @Test
    public void testDiscordantGroups()
    {
        // 3 fragments are required to support a discordant junction, unassigned to other junctions
        List<ReadGroup> discordantCandidates = Lists.newArrayList();

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_2, 5100);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 421, CHR_3, 5000); // unrelated
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 431, CHR_2, 5150);

        List<KnownHotspot> knownHotspots = Lists.newArrayList();
        DiscordantGroups discordantGroups = new DiscordantGroups(
                REGION_1, DEFAULT_MAX_FRAGMENT_LENGTH, knownHotspots, false);

        List<JunctionData> junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);
        assertEquals(0, junctions.size());

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 441, CHR_2, 5050);

        junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);
        assertEquals(2, junctions.size());
        assertTrue(junctions.get(0).junctionGroups().isEmpty());
        assertEquals(540, junctions.get(0).Position);
        assertEquals(FORWARD, junctions.get(0).Orient);

        assertEquals(5050, junctions.get(1).Position);
        assertEquals(REVERSE, junctions.get(1).Orient);

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 540, CHR_2, 105000); // unrelated
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 550, CHR_2, 5300);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 560, CHR_2, 6000); // too far
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 570, CHR_2, 5000);

        junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);
        assertEquals(2, junctions.size());
        assertEquals(8, junctions.get(0).supportingGroups().size());

        // a local DEL needs more support
        discordantCandidates.clear();
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);
        assertEquals(0, junctions.size());

        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        addDiscordantCandidate(discordantCandidates, READ_ID_GENERATOR.nextId(), CHR_1, 401, CHR_1, 1200);
        junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);

        assertEquals(2, junctions.size());
        assertEquals(500, junctions.get(0).Position);
        assertEquals(FORWARD, junctions.get(0).Orient);
        assertEquals(1200, junctions.get(1).Position);
        assertEquals(REVERSE, junctions.get(1).Orient);

    }

    @Test
    public void testDiscordantGroupDistanceFilters()
    {
        // fragments cannot be too far from the junction
        List<ReadGroup> discordantCandidates = Lists.newArrayList();

        PrepRead read1 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2000, CHR_2, 500,true, false, null));
        read1.record().setReadNegativeStrandFlag(true);
        read1.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read1));

        PrepRead read2 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2090, CHR_2, 500,true, false, null));
        read2.record().setReadNegativeStrandFlag(true);
        read2.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read2));

        PrepRead read3 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2180, CHR_2, 500,true, false, null));
        read3.record().setReadNegativeStrandFlag(true);
        read3.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read3));

        // included but its remote region has too little support to form its own junction
        PrepRead read4 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2270, CHR_2, 10000,true, false, null));
        read4.record().setReadNegativeStrandFlag(true);
        read4.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read4));

        PrepRead read5 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2360, CHR_2, 10000,true, false, null));
        read5.record().setReadNegativeStrandFlag(true);
        read5.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read5));

        // ignored for being too far from the junction
        PrepRead read6 = PrepRead.from(createSamRecord(
                READ_ID_GENERATOR.nextId(), CHR_1, 2450, CHR_3, 500,true, false, null));
        read6.record().setReadNegativeStrandFlag(true);
        read6.record().setMateNegativeStrandFlag(false);

        discordantCandidates.add(new ReadGroup(read6));

        List<KnownHotspot> knownHotspots = Lists.newArrayList();

        DiscordantGroups discordantGroups = new DiscordantGroups(
                REGION_1, 500, knownHotspots, false);

        List<JunctionData> junctions = discordantGroups.formDiscordantJunctions(discordantCandidates);
        assertEquals(2, junctions.size());

        assertEquals(5, junctions.get(0).supportingGroups().size());
    }
}
