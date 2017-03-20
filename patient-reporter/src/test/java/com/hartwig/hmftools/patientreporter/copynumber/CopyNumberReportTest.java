package com.hartwig.hmftools.patientreporter.copynumber;

import static org.junit.Assert.assertEquals;

import java.util.List;
import java.util.SortedSet;

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
        final CopyNumberReport first = new CopyNumberReport.Builder().chromosome("18").build();
        final CopyNumberReport second = new CopyNumberReport.Builder().chromosome("Y").build();
        final CopyNumberReport third = new CopyNumberReport.Builder().chromosome("X").build();
        final CopyNumberReport fourth = new CopyNumberReport.Builder().chromosome("12").build();

        final List<CopyNumberReport> reports = Lists.newArrayList(first, second, third, fourth);
        final SortedSet sortedReports = Sets.newTreeSet(reports);

        assertEquals(fourth, sortedReports.first());
        assertEquals(second, sortedReports.last());
    }
}