package com.hartwig.hmftools.viridian.detection;

import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome.MUTUAL;
import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome.NO_CANDIDATES;
import static com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome.ONE_CANDIDATE;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Stream;

import com.hartwig.hmftools.viridian.detection.common.ContigStats;
import com.hartwig.hmftools.viridian.detection.common.SummaryStats;
import com.hartwig.hmftools.viridian.detection.select.ContigRole;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupOutcome;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupRepresentativeSelection;
import com.hartwig.hmftools.viridian.detection.select.OncologyGroupResolution;
import com.hartwig.hmftools.viridian.detection.select.RepresentativeContigCandidate;
import com.hartwig.hmftools.viridian.detection.support.ContigFilterStatus;
import com.hartwig.hmftools.viridian.detection.support.ContigSupport;
import com.hartwig.hmftools.viridian.reference.OncologyGroup;
import com.hartwig.hmftools.viridian.reference.ViralContig;

import org.jetbrains.annotations.Nullable;
import org.junit.Test;

public class DetectedVirusTest
{
    private static final OncologyGroup GROUP_A = new OncologyGroup("Group A");
    private static final OncologyGroup GROUP_B = new OncologyGroup("Group B");
    private static final ViralContig V1 = new ViralContig("v1", 100, "Virus v1", GROUP_A);
    private static final ViralContig V2 = new ViralContig("v2", 100, "Virus v2", GROUP_A);
    private static final ViralContig V3 = new ViralContig("v3", 100, "Virus v3", GROUP_A);
    private static final ViralContig W1 = new ViralContig("w1", 100, "Virus w1", GROUP_B);

    @Test
    public void testFromResolvedGroupCarriesItsRepresentativeMeasurement()
    {
        ContigStats stats = stats(V1);

        List<DetectedVirus> detections = DetectedVirus.from(
                List.of(selection(GROUP_A, ONE_CANDIDATE, V1)), Map.of(V1, stats), Map.of(GROUP_A, 500));

        assertEquals(1, detections.size());
        assertEquals(GROUP_A, detections.get(0).oncologyGroup());
        assertEquals(OncologyGroupResolution.RESOLVED, detections.get(0).resolution());
        assertEquals(stats, detections.get(0).representativeContigStats());
    }

    // The group is detected either way: which of its near-identical strains leads is a separate question, and naming
    // one would override the checks that declined to pick it.
    @Test
    public void testFromUnresolvedGroupIsPresentButNotMeasured()
    {
        List<DetectedVirus> detections = DetectedVirus.from(
                List.of(selection(GROUP_A, MUTUAL, null, List.of(V1, V2), List.of(V3))), Map.of(), Map.of(GROUP_A, 900));

        assertEquals(1, detections.size());
        DetectedVirus detection = detections.get(0);
        assertEquals(OncologyGroupResolution.UNRESOLVED, detection.resolution());
        assertNull(detection.representativeContigStats());

        // Nothing is measured, so these counts are all a reader has to judge the ambiguity by.
        assertEquals(900, detection.groupReadCount());
        assertEquals(3, detection.alignedContigCount());
        assertEquals(2, detection.candidateCount());
        assertEquals(2, detection.comparableCandidateCount());
    }

    @Test
    public void testFromGroupWithNoCandidatesIsAbsent()
    {
        List<DetectedVirus> detections =
                DetectedVirus.from(List.of(selection(GROUP_B, NO_CANDIDATES, null)), Map.of(), Map.of());

        assertTrue(detections.isEmpty());
    }

    // A representative with no measurement means the stages disagree about what was assigned.
    @Test
    public void testFromRepresentativeMissingItsMeasurementRejected()
    {
        List<OncologyGroupRepresentativeSelection> selections = List.of(selection(GROUP_A, ONE_CANDIDATE, V1));

        assertThrows(
                IllegalStateException.class,
                () -> DetectedVirus.from(selections, Map.of(W1, stats(W1)), Map.of(GROUP_A, 500)));
    }

    private static ContigStats stats(ViralContig contig)
    {
        return new ContigStats(
                contig, 10, 1, 80, SummaryStats.from(new int[contig.length()]), SummaryStats.from(new int[] { 100 }));
    }

    private static OncologyGroupRepresentativeSelection selection(
            OncologyGroup group, OncologyGroupOutcome outcome, @Nullable ViralContig representative, ViralContig... others)
    {
        return selection(group, outcome, representative, List.of(others), List.of());
    }

    private static OncologyGroupRepresentativeSelection selection(
            OncologyGroup group, OncologyGroupOutcome outcome, @Nullable ViralContig representative,
            List<ViralContig> others, List<ViralContig> rejected)
    {
        List<RepresentativeContigCandidate> candidates = Stream.concat(
                        representative != null ? Stream.of(representative) : Stream.empty(), others.stream())
                .map(contig -> candidate(contig, contig.equals(representative)))
                .toList();
        List<ContigSupport> rejectedSupports = rejected.stream()
                .map(contig -> new ContigSupport(stats(contig), ContigFilterStatus.LOW_COVERAGE, 0, null, 0.1))
                .toList();

        return new OncologyGroupRepresentativeSelection(group, outcome, candidates, rejectedSupports);
    }

    private static RepresentativeContigCandidate candidate(ViralContig contig, boolean representative)
    {
        ContigSupport support = new ContigSupport(stats(contig), ContigFilterStatus.CANDIDATE, 0, null, 1.0);
        return new RepresentativeContigCandidate(
                support, 1, true, Set.of(), Set.of(),
                representative ? ContigRole.REPRESENTATIVE : ContigRole.SECONDARY);
    }
}
