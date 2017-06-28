package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import org.junit.Test;

public class CopyNumberReportTest {

    @Test
    public void canResolveType() {
        final CopyNumberReport gain = new CopyNumberReport.Builder().copyNumber(5).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_GAIN, gain.resolveType());

        final CopyNumberReport loss = new CopyNumberReport.Builder().copyNumber(0).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_LOSS, loss.resolveType());

        final CopyNumberReport neutral = new CopyNumberReport.Builder().copyNumber(2).build();
        assertEquals(CopyNumberReport.COPY_NUMBER_NEUTRAL, neutral.resolveType());
    }

    @Test
    public void canSort() {
        final CopyNumberReport first = new CopyNumberReport.Builder().chromosome("18").chromosomeBand("p2").build();
        final CopyNumberReport second = new CopyNumberReport.Builder().chromosome("Y").chromosomeBand("q1").build();
        final CopyNumberReport third = new CopyNumberReport.Builder().chromosome("X").chromosomeBand("q1").build();
        final CopyNumberReport fourth = new CopyNumberReport.Builder().chromosome("12").chromosomeBand("q1").build();
        final CopyNumberReport fifth = new CopyNumberReport.Builder().chromosome("18").chromosomeBand("p1").build();
        final CopyNumberReport sixth = new CopyNumberReport.Builder().chromosome("X").chromosomeBand("p1").build();

        final List<CopyNumberReport> reports = Lists.newArrayList(first, second, third, fourth, fifth, sixth);
        final List<CopyNumberReport> sortedReports = Lists.newArrayList(Sets.newTreeSet(reports));

        assertEquals(fourth, sortedReports.get(0));
        assertEquals(fifth, sortedReports.get(1));
        assertEquals(first, sortedReports.get(2));
        assertEquals(sixth, sortedReports.get(3));
        assertEquals(third, sortedReports.get(4));
        assertEquals(second, sortedReports.get(5));
    }
}