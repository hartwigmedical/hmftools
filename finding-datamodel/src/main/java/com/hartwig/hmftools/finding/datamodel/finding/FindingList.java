package com.hartwig.hmftools.finding.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.finding.datamodel.IFindingList;
import com.hartwig.hmftools.finding.datamodel.RecordBuilder;

import jakarta.validation.constraints.NotNull;

@RecordBuilder
public record FindingList<T extends Finding>(
        @NotNull FindingsStatus status,
        @NotNull List<T> findings) implements IFindingList<T>
{

}
