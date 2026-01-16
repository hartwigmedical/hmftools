package com.hartwig.hmftools.datamodel.finding;

import java.util.List;
import java.util.Set;

import com.hartwig.hmftools.datamodel.driver.DriverSource;

public class FindingsQuery<T extends Driver>
{

    private final List<T> findings;
    private Set<DriverSource> driverSources = null;
    private Set<ReportedStatus> reportedStatuses = null;

    public FindingsQuery(List<T> findings)
    {
        this.findings = findings;
    }

    public FindingsQuery<T> driverSources(DriverSource... driverSources)
    {
        this.driverSources = Set.of(driverSources);
        return this;
    }

    public FindingsQuery<T> germlineOnly()
    {
        return driverSources(DriverSource.GERMLINE);
    }

    public FindingsQuery<T> somaticOnly()
    {
        return driverSources(DriverSource.SOMATIC);
    }

    public FindingsQuery<T> reportedStatuses(ReportedStatus... reportedStatuses)
    {
        this.reportedStatuses = Set.of(reportedStatuses);
        return this;
    }

    public List<T> results()
    {
        return findings.stream()
                .filter(o -> driverSources == null || driverSources.contains(o.driverSource()))
                .filter(o -> reportedStatuses == null || reportedStatuses.contains(o.reportedStatus())).toList();

    }
}
