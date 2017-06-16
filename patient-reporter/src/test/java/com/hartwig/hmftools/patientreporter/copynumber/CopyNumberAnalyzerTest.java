package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.Set;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.copynumber.CopyNumber;
import com.hartwig.hmftools.common.copynumber.freec.ImmutableFreecCopyNumber;
import com.hartwig.hmftools.common.region.hmfslicer.HmfGenomeRegion;
import com.hartwig.hmftools.common.region.hmfslicer.ImmutableHmfGenomeRegion;

import org.jetbrains.annotations.NotNull;
import org.junit.Test;

public class CopyNumberAnalyzerTest {

    private static final String CHROMOSOME = "X";
    private static final String ENTREZ_ID = "11";
    private static final int TRANSCRIPT_VERSION = 1;
    private static final double EPSILON = 1.0e-10;

    @Test
    public void worksAsExpected() {
        final HmfGenomeRegion first = region(101, 200, "TRANS1", "GENE1", "p1");
        final HmfGenomeRegion second = region(301, 400, "TRANS2", "GENE2", "p2");
        final HmfGenomeRegion third = region(601, 700, "TRANS3", "GENE3", "p3");

        final Set<HmfGenomeRegion> regions = Sets.newHashSet(first, second, third);

        final List<CopyNumber> copyNumbers = Lists.newArrayList();
        copyNumbers.add(copyNumber(11, 20, 1));
        copyNumbers.add(copyNumber(121, 170, 4));
        copyNumbers.add(copyNumber(191, 260, 5));
        copyNumbers.add(copyNumber(291, 500, 0));

        final CopyNumberAnalyzer copyNumberAnalyzer = new CopyNumberAnalyzer(regions);
        final CopyNumberAnalysis analysis = copyNumberAnalyzer.run(copyNumbers);

        final List<CopyNumberReport> findings = analysis.findings();
        assertEquals(1, findings.size());
        assertEquals(CHROMOSOME, findings.get(0).chromosome());
        assertEquals(second.gene(), findings.get(0).gene());
        assertEquals(second.transcript(), findings.get(0).transcript());
        assertEquals(second.chromosomeBand(), findings.get(0).chromosomeBand());
        assertEquals(0, findings.get(0).copyNumber());

        final CopyNumberStats firstStat = analysis.stats().get(first);
        assertEquals(2, firstStat.min());
        assertEquals(5, firstStat.max());
        assertEquals(3.3, firstStat.mean(), EPSILON);

        final CopyNumberStats secondStat = analysis.stats().get(second);
        assertEquals(0, secondStat.min());
        assertEquals(0, secondStat.max());
        assertEquals(0D, secondStat.mean(), EPSILON);

        final CopyNumberStats thirdStat = analysis.stats().get(third);
        assertEquals(2, thirdStat.min());
        assertEquals(2, thirdStat.max());
        assertEquals(2D, thirdStat.mean(), EPSILON);
    }

    @NotNull
    private static HmfGenomeRegion region(final long start, final long end, @NotNull final String transcriptId,
            @NotNull final String gene, @NotNull final String chromosomeBand) {
        return new ImmutableHmfGenomeRegion(CHROMOSOME, start, end, transcriptId, TRANSCRIPT_VERSION, gene,
                chromosomeBand, ENTREZ_ID);
    }

    @NotNull
    private static CopyNumber copyNumber(final long start, final long end, final int value) {
        return ImmutableFreecCopyNumber.builder().chromosome(CHROMOSOME).start(start).end(end).value(value).build();
    }
}
