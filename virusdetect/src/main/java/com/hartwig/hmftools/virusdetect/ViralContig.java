package com.hartwig.hmftools.virusdetect;

// TODO: maybe rename to "viral genome"? "contig" is not that descriptive. But that's a big change
// A contig which is a virus genome.
public record ViralContig(
        // Original contig name in the viral reference.
        String name,
        int length,
        // Human-readable name of the virus. Only for readability purposes.
        String virusName,
        OncologyGroup oncologyGroup
)
{
}
