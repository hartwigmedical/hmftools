package com.hartwig.hmftools.serve.extraction.hotspot.tools;

import java.util.Set;

import com.google.common.collect.Sets;

import org.jetbrains.annotations.NotNull;

final class AnnotatedHotspotCurationFactory {

    public static final Set<String> RETIRED_TRANSCRIPTS = createRetiredTranscriptWhitelist();
    public static final Set<String> CHANGED_TRANSCRIPTS = createChangedTranscriptWhitelist();

    @NotNull
    private static Set<String> createRetiredTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This transcript on RHOA is coding in v37 but non-coding in v38
        transcripts.add("ENST00000431929");

        // This transcript on CDKN2A has retired in v38
        transcripts.add("ENST00000361570");

        // This transcript on BCL2L12 has completely changed in 38 and can't match anymore.
        transcripts.add("ENST00000246784");

        return transcripts;
    }

    @NotNull
    private static Set<String> createChangedTranscriptWhitelist() {
        Set<String> transcripts = Sets.newHashSet();

        // This MYC transcript is v2 in 37 and v6 in 38 and they differ by 15 AA.
        transcripts.add("ENST00000377970");

        // This TCF7L2 transcript in v1 in 37 and v5 in 38 and they differ by 5 AA.
        transcripts.add("ENST00000543371");

        // This FGFR2 transcript is v6 in 37 and v10 in 38 and they differ by 36 AA (spread across transcript).
        transcripts.add("ENST00000351936");

        // This MED12 transcript is v6 in 37 and v10 in 38 and they differ with 53 AA (spread across transcript)
        transcripts.add("ENST00000333646");

        return transcripts;
    }

    private AnnotatedHotspotCurationFactory() {
    }
}
