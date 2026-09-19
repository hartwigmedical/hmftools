package com.hartwig.hmftools.virusdetect;

// A contig the prefilter rejected, with the criterion it failed. Reported only: nothing acts on the reason.
public record RejectedContig(
        ContigSupport stats,
        ContigFilterStatus reason
)
{
}
