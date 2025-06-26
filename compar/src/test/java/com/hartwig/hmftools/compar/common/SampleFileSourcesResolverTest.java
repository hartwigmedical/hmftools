package com.hartwig.hmftools.compar.common;

import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.checkAddDirSeparator;

import static org.junit.Assert.assertEquals;

import com.hartwig.pipeline.reference.api.PipelineOutputStructure;
import com.hartwig.pipeline.reference.api.PipelineVersion;

import org.junit.Test;

public class SampleFileSourcesResolverTest
{
    static final String WILDCARD_SAMPLE_DIR = "/base/dir/*/";
    static final String TUMOR_SAMPLE_ID = "tumor123";
    static final String REF_SAMPLE_ID = "ref123";
    static final String SAMPLE_DIR = "/base/dir/tumor123/";

    @Test
    public void testOnlySampleDirArgument()
    {
        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, "", "", "", null, null);
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertPipeline5Structure(victim, SAMPLE_DIR, "", "");
    }

    @Test
    public void testToolAndSampleDirArguments()
    {
        String chordPath = "chord/path/*/";
        String unfilteredSomaticVcfPath = "/path/to/*/unfiltered";
        String somaticVcfPath = "/path/to/$/somatic_vcf";

        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, chordPath, somaticVcfPath, unfilteredSomaticVcfPath, null, null);
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertEquals(SAMPLE_DIR + "chord/path/" + TUMOR_SAMPLE_ID + "/", victim.chord());
        assertEquals("/path/to/" + TUMOR_SAMPLE_ID + "/unfiltered", victim.somaticUnfilteredVcf());
        assertEquals("/path/to/" + REF_SAMPLE_ID + "/somatic_vcf", victim.somaticVcf());
    }

    @Test
    public void testToolOnlyArguments()
    {
        String chordPath = "chord/path/*/";
        String unfilteredSomaticVcfPath = "/path/to/*/unfiltered";
        String somaticVcfPath = "/path/to/$/somatic_vcf";

        FileSources fileSources = createVictim("", chordPath, somaticVcfPath, unfilteredSomaticVcfPath, null, null);
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertEquals("chord/path/" + TUMOR_SAMPLE_ID + "/", victim.chord());
        assertEquals("/path/to/" + TUMOR_SAMPLE_ID + "/unfiltered", victim.somaticUnfilteredVcf());
        assertEquals("/path/to/" + REF_SAMPLE_ID + "/somatic_vcf", victim.somaticVcf());
    }

    @Test
    public void testPipeline5V60Format()
    {
        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, "", "", "", "6.0", "PIPELINE5");
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertPipeline5Structure(victim, SAMPLE_DIR, "", SAMPLE_DIR + "sage_somatic/" + TUMOR_SAMPLE_ID + ".sage.somatic.vcf.gz");
    }

    @Test
    public void testPipeline5V534Format()
    {
        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, "", "", "", "5.34", "PIPELINE5");
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertPipeline5Structure(victim, SAMPLE_DIR, "", SAMPLE_DIR + "sage_somatic/" + TUMOR_SAMPLE_ID + ".sage.somatic.vcf.gz");
    }

    @Test
    public void testOncoAnalyserV2Format()
    {
        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, "", "", "", "6.0", "ONCOANALYSER");
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertOncoAnalyserStructure(victim, SAMPLE_DIR, "", SAMPLE_DIR + "sage/somatic/" + TUMOR_SAMPLE_ID + ".sage.somatic.vcf.gz");
    }

    @Test
    public void testDatabaseV60Format()
    {
        FileSources fileSources = createVictim(WILDCARD_SAMPLE_DIR, "", "", "", "6.0", "DATABASE");
        SampleFileSources victim = SampleFileSources.fromFileSources(fileSources, TUMOR_SAMPLE_ID, REF_SAMPLE_ID);

        assertDatabaseStructure(victim, SAMPLE_DIR, "", SAMPLE_DIR + "sage/somatic/" + TUMOR_SAMPLE_ID + ".sage.somatic.vcf.gz");
    }

    private static FileSources createVictim(String sampleDir, String chordDir, String somaticVcf, String unfilteredSomaticVcf,
            String pipelineVersion, String pipelineOutputStructure)
    {
        return new FileSources("ref", sampleDir, "", "", "", "", "", chordDir,
                "", "", somaticVcf, unfilteredSomaticVcf, "", "", "",
                "", "", "", "",
                pipelineVersion != null ? PipelineVersion.fromString(pipelineVersion) : null,
                pipelineOutputStructure != null ? PipelineOutputStructure.valueOf(pipelineOutputStructure) : null);
    }

    private static void assertPipeline5Structure(final SampleFileSources victim, final String sampleDir, final String somaticVcf,
            final String somaticUnfilteredVcf)
    {
        String sampleDirWithTrailingSlash = checkAddDirSeparator(sampleDir);
        assertEquals(sampleDirWithTrailingSlash + "linx/", victim.linx());
        assertEquals(sampleDirWithTrailingSlash + "linx_germline/", victim.linxGermline());
        assertEquals(sampleDirWithTrailingSlash + "purple/", victim.purple());
        assertEquals(sampleDirWithTrailingSlash + "cuppa/", victim.cuppa());
        assertEquals(sampleDirWithTrailingSlash + "lilac/", victim.lilac());
        assertEquals(sampleDirWithTrailingSlash + "chord/", victim.chord());
        assertEquals(sampleDirWithTrailingSlash + "peach/", victim.peach());
        assertEquals(sampleDirWithTrailingSlash + "virusintrprtr/", victim.virus());
        assertEquals(sampleDirWithTrailingSlash + TUMOR_SAMPLE_ID + "/bam_metrics/", victim.tumorFlagstat());
        assertEquals(sampleDirWithTrailingSlash + REF_SAMPLE_ID + "/bam_metrics/", victim.germlineFlagstat());
        assertEquals(sampleDirWithTrailingSlash + TUMOR_SAMPLE_ID + "/bam_metrics/", victim.tumorBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + REF_SAMPLE_ID + "/bam_metrics/", victim.germlineBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + REF_SAMPLE_ID + "/snp_genotype/", victim.snpGenotype());
        assertEquals(sampleDirWithTrailingSlash + "cider/", victim.cider());
        assertEquals(sampleDirWithTrailingSlash + "teal/", victim.teal());

        assertEquals(somaticVcf, victim.somaticVcf());
        assertEquals(somaticUnfilteredVcf, victim.somaticUnfilteredVcf());
    }

    private static void assertOncoAnalyserStructure(final SampleFileSources victim, final String sampleDir, final String somaticVcf,
            final String somaticUnfilteredVcf)
    {
        String sampleDirWithTrailingSlash = checkAddDirSeparator(sampleDir);
        assertEquals(sampleDirWithTrailingSlash + "linx/somatic_annotations/", victim.linx());
        assertEquals(sampleDirWithTrailingSlash + "linx/germline_annotations/", victim.linxGermline());
        assertEquals(sampleDirWithTrailingSlash + "purple/", victim.purple());
        assertEquals(sampleDirWithTrailingSlash + "cuppa/", victim.cuppa());
        assertEquals(sampleDirWithTrailingSlash + "lilac/", victim.lilac());
        assertEquals(sampleDirWithTrailingSlash + "chord/", victim.chord());
        assertEquals(sampleDirWithTrailingSlash + "peach/", victim.peach());
        assertEquals(sampleDirWithTrailingSlash + "virusinterpreter/", victim.virus());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + TUMOR_SAMPLE_ID + "_bamtools/", victim.tumorFlagstat());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + REF_SAMPLE_ID + "_bamtools/", victim.germlineFlagstat());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + TUMOR_SAMPLE_ID + "_bamtools/", victim.tumorBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + REF_SAMPLE_ID + "_bamtools/", victim.germlineBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + "cider/", victim.cider());
        assertEquals(sampleDirWithTrailingSlash + "teal/", victim.teal());

        // Not available in OA, so default to "base default"
        assertEquals(sampleDirWithTrailingSlash + REF_SAMPLE_ID + "/snp_genotype/", victim.snpGenotype());

        assertEquals(somaticVcf, victim.somaticVcf());
        assertEquals(somaticUnfilteredVcf, victim.somaticUnfilteredVcf());
    }

    private static void assertDatabaseStructure(final SampleFileSources victim, final String sampleDir, final String somaticVcf,
            final String somaticUnfilteredVcf)
    {
        String sampleDirWithTrailingSlash = checkAddDirSeparator(sampleDir);
        assertEquals(sampleDirWithTrailingSlash + "linx/somatic_annotations/", victim.linx());
        assertEquals(sampleDirWithTrailingSlash + "linx/germline_annotations/", victim.linxGermline());
        assertEquals(sampleDirWithTrailingSlash + "purple/", victim.purple());
        assertEquals(sampleDirWithTrailingSlash + "cuppa/", victim.cuppa());
        assertEquals(sampleDirWithTrailingSlash + "lilac/", victim.lilac());
        assertEquals(sampleDirWithTrailingSlash + "chord/", victim.chord());
        assertEquals(sampleDirWithTrailingSlash + "peach/", victim.peach());
        assertEquals(sampleDirWithTrailingSlash + "virusinterpreter/", victim.virus());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + TUMOR_SAMPLE_ID + "_bamtools/", victim.tumorFlagstat());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + REF_SAMPLE_ID + "_bamtools/", victim.germlineFlagstat());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + TUMOR_SAMPLE_ID + "_bamtools/", victim.tumorBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + "bamtools/" + REF_SAMPLE_ID + "_bamtools/", victim.germlineBamMetrics());
        assertEquals(sampleDirWithTrailingSlash + "snp_genotype/" + REF_SAMPLE_ID + "/", victim.snpGenotype());
        assertEquals(sampleDirWithTrailingSlash + "cider/", victim.cider());
        assertEquals(sampleDirWithTrailingSlash + "teal/", victim.teal());

        assertEquals(somaticVcf, victim.somaticVcf());
        assertEquals(somaticUnfilteredVcf, victim.somaticUnfilteredVcf());
    }
}
