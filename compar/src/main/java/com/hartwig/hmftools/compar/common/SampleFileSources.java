package com.hartwig.hmftools.compar.common;

public record SampleFileSources(String source, String linx, String purple, String linxGermline, String cuppa, String lilac, String chord,
                                String peach, String virus, String somaticVcf, String somaticUnfilteredVcf, String tumorFlagstat,
                                String germlineFlagstat, String tumorBamMetrics, String germlineBamMetrics, String snpGenotype,
                                String cider, String teal)
{
    public static SampleFileSources fromFileSources(final FileSources fileSources, final String sampleId, final String germlineSampleId)
    {
        SampleFileSourceResolver resolver = new SampleFileSourceResolver(fileSources, sampleId, germlineSampleId);

        return new SampleFileSources(
                fileSources.source(),
                resolver.resolveLinxDirectory(),
                resolver.resolvePurpleDirectory(),
                resolver.resolveLinxGermlineDirectory(),
                resolver.resolveCuppaDirectory(),
                resolver.resolveLilacDirectory(),
                resolver.resolveChordDirectory(),
                resolver.resolvePeachDirectory(),
                resolver.resolveVirusDirectory(),
                resolver.resolveSomaticVcfPath(),
                resolver.resolveSomaticUnfilteredVcfPath(),
                resolver.resolveTumorFlagstatDirectory(),
                resolver.resolveGermlineFlagstatDirectory(),
                resolver.resolveTumorBamMetricsDirectory(),
                resolver.resolveGermlineBamMetricsDirectory(),
                resolver.resolveSnpGenotypeDirectory(),
                resolver.resolveCiderDirectory(),
                resolver.resolveTealDirectory());
    }
}
