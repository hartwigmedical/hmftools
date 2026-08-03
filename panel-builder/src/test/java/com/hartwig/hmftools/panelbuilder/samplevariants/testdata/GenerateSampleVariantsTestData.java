package com.hartwig.hmftools.panelbuilder.samplevariants.testdata;

import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.REF_GENOME;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.addRefGenomeFile;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_FWD;
import static com.hartwig.hmftools.common.genome.region.Orientation.ORIENT_REV;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE_DESC;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.OUTPUT_DIR_DESC;

import java.util.List;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.purple.GermlineStatus;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.common.variant.CodingEffect;
import com.hartwig.hmftools.panelbuilder.samplevariants.testdata.SampleVariantsTestData.GermlineSv;
import com.hartwig.hmftools.panelbuilder.samplevariants.testdata.SampleVariantsTestData.SmallVariant;
import com.hartwig.hmftools.panelbuilder.samplevariants.testdata.SampleVariantsTestData.Sv;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

// Command line tool that writes a curated set of synthetic sample variant files against a real reference genome, for running
// PanelBuilder's sample variant feature. See SAMPLE_VARIANTS_TEST_DATA.md. No patient data is used.
public class GenerateSampleVariantsTestData
{
    private static final String APP_NAME = "GenerateSampleVariantsTestData";
    private static final String DEFAULT_SAMPLE = "FAKE01T";

    private static final Logger LOGGER = LogManager.getLogger(GenerateSampleVariantsTestData.class);

    public static void main(final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(APP_NAME);
        addRefGenomeFile(configBuilder, true);
        configBuilder.addConfigItem(OUTPUT_DIR, true, OUTPUT_DIR_DESC);
        configBuilder.addConfigItem(SAMPLE, false, SAMPLE_DESC);
        addLoggingOptions(configBuilder);
        if(!configBuilder.parseCommandLine(args))
        {
            System.exit(1);
        }
        setLogLevel(configBuilder);

        String refGenomePath = configBuilder.getValue(REF_GENOME);
        String outputDir = configBuilder.getValue(OUTPUT_DIR);
        String sampleId = configBuilder.getValue(SAMPLE, DEFAULT_SAMPLE);

        RefGenomeInterface genome = RefGenomeSource.loadRefGenome(refGenomePath);

        SampleVariantsTestData.OutputDirs dirs = SampleVariantsTestData.generate(genome, sampleId, outputDir, curatedGrch38Specs());
        LOGGER.info("Wrote sample variant test data for {}", sampleId);
        LOGGER.info("  purple_dir        {}", dirs.purpleDir());
        LOGGER.info("  linx_dir          {}", dirs.linxDir());
        LOGGER.info("  linx_germline_dir {}", dirs.linxGermlineDir());
    }

    // Public cancer variants placed at approximate GRCh38 coordinates. One entry per PanelBuilder sample variant probe type.
    private static SampleVariantsTestData.Specs curatedGrch38Specs()
    {
        List<SmallVariant> somaticSmall = List.of(
                // Somatic SNV driver (BRAF V600E locus).
                SmallVariant.snv("chr7", 140753336, true, "BRAF", CodingEffect.MISSENSE, GermlineStatus.DIPLOID, 0.0, 70, 30),
                // Somatic SNV non-driver (PIK3CA locus), passes the non-driver filters.
                SmallVariant.snv("chr3", 179234297, false, "PIK3CA", CodingEffect.MISSENSE, GermlineStatus.DIPLOID, 0.0, 70, 30),
                // Somatic small deletion non-driver (TP53 locus).
                SmallVariant.deletion("chr17", 7675000, 3, false, "TP53", CodingEffect.NONE, GermlineStatus.DIPLOID, 70, 30));

        List<SmallVariant> germlineSmall = List.of(
                // Germline SNV driver (BRCA2 locus).
                SmallVariant.snv("chr13", 32340000, true, "BRCA2", CodingEffect.MISSENSE, GermlineStatus.HET_DELETION, 0.0, 60, 40));

        List<Sv> somaticSvs = List.of(
                // SV fusion driver (TMPRSS2-ERG, intra-chromosomal on chr21).
                new Sv(
                        "chr21", 41498119, ORIENT_FWD, "chr21", 38445621, ORIENT_REV, "",
                        Sv.DriverKind.FUSION, "TMPRSS2_ERG", 0.3, 1.0, 1.0, 30),
                // SV amplification driver (ERBB2 locus).
                new Sv(
                        "chr17", 39700000, ORIENT_REV, "chr17", 39730000, ORIENT_FWD, "",
                        Sv.DriverKind.AMP, "ERBB2", 0.3, 3.0, 3.0, 30),
                // SV deletion driver (CDKN2A locus).
                new Sv(
                        "chr9", 21967000, ORIENT_FWD, "chr9", 21995000, ORIENT_REV, "",
                        Sv.DriverKind.DEL, "CDKN2A", 0.3, 1.0, 1.0, 30),
                // SV disruption driver. Breakends inside the large BRCA2 exon 11 (high mappability); PTEN was avoided as its PTENP1
                // pseudogene makes probes there low quality.
                new Sv(
                        "chr13", 32333000, ORIENT_FWD, "chr13", 32338000, ORIENT_REV, "",
                        Sv.DriverKind.DISRUPTION, "BRCA2", 0.3, 1.0, 1.0, 20));

        List<GermlineSv> germlineSvs = List.of(
                // Germline SV driver (BRCA1 locus).
                new GermlineSv("chr17", 43044000, ORIENT_FWD, "chr17", 43120000, ORIENT_REV, "", "BRCA1"));

        return new SampleVariantsTestData.Specs(somaticSmall, germlineSmall, somaticSvs, germlineSvs);
    }
}
