package com.hartwig.hmftools.datamodel.finding;

import java.util.List;

import com.hartwig.hmftools.datamodel.purple.PurpleMicrosatelliteStatus;

import org.jetbrains.annotations.NotNull;

public record MicrosatelliteStability(
        @NotNull String findingKey,
        double microsatelliteIndelsPerMb,
        @NotNull PurpleMicrosatelliteStatus microsatelliteStatus,
        @NotNull List<GainDeletion> lohCopyNumbers,
        @NotNull List<String> genes)
implements Finding {}
