package com.hartwig.hmftools.finding.util;

import java.util.Comparator;
import java.util.List;
import java.util.Objects;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import com.hartwig.hmftools.finding.datamodel.FindingRecord;
import com.hartwig.hmftools.finding.datamodel.driver.Driver;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingList;
import com.hartwig.hmftools.finding.datamodel.driver.DriverFindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.Finding;
import com.hartwig.hmftools.finding.datamodel.finding.FindingList;
import com.hartwig.hmftools.finding.datamodel.finding.FindingListBuilder;
import com.hartwig.hmftools.finding.datamodel.finding.FindingStatus;
import com.hartwig.hmftools.finding.datamodel.finding.IFindingList;

import org.jspecify.annotations.Nullable;

import jakarta.validation.constraints.NotNull;

public class FindingRecordConverterUtil
{
    public static Function<FindingRecord, FindingRecord> listConverter(List<Function<FindingRecord, FindingRecord>> converters)
    {
        return record -> converters.stream().reduce(record, (r, c) -> c.apply(r), (r1, r2) -> r1);
    }

    @NotNull
    public static <I extends Finding, O extends Finding> FindingList<O> convertFindingList(@NotNull IFindingList<I> findingList,
            Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @Nullable Function<I, O> findingConverter,
            @Nullable Comparator<O> comparator)
    {
        return FindingListBuilder.<O>builder()
                .status(findingsStatusConverter.apply(findingList.status()))
                .findings(convert(findingList.findings(), findingConverter, comparator))
                .build();
    }

    @NotNull
    public static <I extends Driver, O extends Driver> DriverFindingList<O> convertDriverFindingList(
            @NotNull IFindingList<I> driverFindingList,
            Function<FindingStatus, FindingStatus> findingsStatusConverter,
            @Nullable Function<I, O> findingConverter,
            @Nullable Comparator<O> comparator)
    {
        return DriverFindingListBuilder.<O>builder()
                .status(findingsStatusConverter.apply(driverFindingList.status()))
                .findings(convert(driverFindingList.findings(), findingConverter, comparator))
                .build();
    }

    static <I, O> List<O> convert(List<I> list, Function<I, O> converter, @Nullable Comparator<O> comparator)
    {
        Stream<O> stream = list.stream().map(converter).filter(Objects::nonNull);
        if(comparator != null)
        {
            stream = stream.sorted(comparator);
        }
        return stream.collect(Collectors.toList());
    }
}
