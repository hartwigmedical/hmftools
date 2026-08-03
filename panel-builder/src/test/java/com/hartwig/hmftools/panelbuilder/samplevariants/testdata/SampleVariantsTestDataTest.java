package com.hartwig.hmftools.panelbuilder.samplevariants.testdata;

import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.test.MockRefGenome;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.panelbuilder.TargetMetadata;
import com.hartwig.hmftools.panelbuilder.samplevariants.GermlineMutation;
import com.hartwig.hmftools.panelbuilder.samplevariants.GermlineSv;
import com.hartwig.hmftools.panelbuilder.samplevariants.SomaticMutation;
import com.hartwig.hmftools.panelbuilder.samplevariants.SomaticSv;
import com.hartwig.hmftools.panelbuilder.samplevariants.testdata.SampleVariantsTestData.Sv;

import org.junit.Test;

public class SampleVariantsTestDataTest
{
    private static final String SAMPLE_ID = "FAKE01T";
    private static final List<String> CHROMOSOMES = List.of("chr7", "chr9", "chr10", "chr17", "chr21");
    private static final int CHROMOSOME_LENGTH = 20000;

    private static MockRefGenome mockGenome()
    {
        MockRefGenome genome = new MockRefGenome(true);
        genome.ChromosomeLengths.clear();
        genome.RefGenomeMap.clear();
        for(String chromosome : CHROMOSOMES)
        {
            genome.ChromosomeLengths.put(chromosome, CHROMOSOME_LENGTH);
            genome.RefGenomeMap.put(chromosome, MockRefGenome.generateRandomBases(CHROMOSOME_LENGTH));
        }
        return genome;
    }

    private static SampleVariantsTestData.Specs specs()
    {
        List<SampleVariantsTestData.SmallVariant> somaticSmall = List.of(
                SampleVariantsTestData.SmallVariant.snv(
                        "chr7", 5000, true, "BRAF", CodingEffect.MISSENSE, GermlineStatus.DIPLOID, 0.0, 70, 30),
                SampleVariantsTestData.SmallVariant.snv(
                        "chr9", 5000, false, "", CodingEffect.NONE, GermlineStatus.DIPLOID, 0.0, 70, 30),
                SampleVariantsTestData.SmallVariant.deletion(
                        "chr9", 6000, 3, false, "", CodingEffect.NONE, GermlineStatus.DIPLOID, 70, 30));

        List<SampleVariantsTestData.SmallVariant> germlineSmall = List.of(
                SampleVariantsTestData.SmallVariant.snv(
                        "chr17", 5000, true, "BRCA1", CodingEffect.MISSENSE, GermlineStatus.HET_DELETION, 0.0, 60, 40));

        List<Sv> somaticSvs = List.of(
                new Sv(
                        "chr21", 5000, ORIENT_FWD, "chr7", 8000, ORIENT_REV, "",
                        Sv.DriverKind.FUSION, "TMPRSS2_ERG", 0.3, 1.0, 1.0, 30),
                new Sv(
                        "chr17", 7000, ORIENT_REV, "chr17", 9000, ORIENT_FWD, "",
                        Sv.DriverKind.AMP, "ERBB2", 0.3, 1.5, 1.5, 30),
                new Sv(
                        "chr9", 8000, ORIENT_FWD, "chr9", 9000, ORIENT_REV, "",
                        Sv.DriverKind.DEL, "CDKN2A", 0.3, 1.0, 1.0, 30),
                new Sv(
                        "chr10", 5000, ORIENT_FWD, "chr10", 6000, ORIENT_REV, "",
                        Sv.DriverKind.DISRUPTION, "PTEN", 0.3, 1.0, 1.0, 20));

        List<SampleVariantsTestData.GermlineSv> germlineSvs = List.of(
                new SampleVariantsTestData.GermlineSv("chr17", 8000, ORIENT_FWD, "chr17", 8500, ORIENT_REV, "", "BRCA1"));

        return new SampleVariantsTestData.Specs(somaticSmall, germlineSmall, somaticSvs, germlineSvs);
    }

    @Test
    public void generatesLoadableSampleVariants() throws IOException
    {
        String outputRoot = Files.createTempDirectory("sample_variants").toString();
        SampleVariantsTestData.OutputDirs dirs = SampleVariantsTestData.generate(mockGenome(), SAMPLE_ID, outputRoot, specs());

        List<SomaticMutation> somaticMutations = SomaticMutation.load(SAMPLE_ID, dirs.purpleDir());
        assertEquals(3, somaticMutations.size());
        assertEquals(1, somaticMutations.stream().filter(SomaticMutation::isDriver).count());

        List<GermlineMutation> germlineMutations = GermlineMutation.load(SAMPLE_ID, dirs.purpleDir());
        assertEquals(1, germlineMutations.size());
        assertEquals(TargetMetadata.Type.SAMPLE_GERMLINE_SNV_INDEL_DRIVER, germlineMutations.get(0).targetType());

        List<SomaticSv> somaticSvs = SomaticSv.load(SAMPLE_ID, dirs.purpleDir(), dirs.linxDir());
        Set<TargetMetadata.Type> svTypes = somaticSvs.stream().map(SomaticSv::targetType).collect(Collectors.toSet());
        assertEquals(4, somaticSvs.size());
        assertTrue(svTypes.contains(TargetMetadata.Type.SAMPLE_SV_FUSION_DRIVER));
        assertTrue(svTypes.contains(TargetMetadata.Type.SAMPLE_SV_AMP_DRIVER));
        assertTrue(svTypes.contains(TargetMetadata.Type.SAMPLE_SV_DEL_DRIVER));
        assertTrue(svTypes.contains(TargetMetadata.Type.SAMPLE_SV_DISRUPTION_DRIVER));

        List<GermlineSv> germlineSvs = GermlineSv.load(SAMPLE_ID, dirs.linxGermlineDir());
        assertEquals(1, germlineSvs.size());
        assertEquals(TargetMetadata.Type.SAMPLE_GERMLINE_SV_DRIVER, germlineSvs.get(0).targetType());

        // Building the probe definitions must not throw for any loaded variant.
        somaticMutations.forEach(SomaticMutation::generateProbe);
        germlineMutations.forEach(GermlineMutation::generateProbe);
        somaticSvs.forEach(SomaticSv::generateProbe);
        germlineSvs.forEach(GermlineSv::generateProbe);
    }
}
