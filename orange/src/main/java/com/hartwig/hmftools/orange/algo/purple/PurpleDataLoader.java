package com.hartwig.hmftools.orange.algo.purple;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.Nullable;

public final class PurpleDataLoader
{
    public static PurpleData load(final OrangeConfig config, final Map<String,DriverGene> driverGenes) throws IOException
    {
        String tumorSample = config.tumorSampleId();
        String referenceSample = config.wgsRefConfig() != null ? config.wgsRefConfig().referenceSampleId() : null;
        String rnaSample = config.rnaConfig() != null ? config.rnaConfig().rnaSampleId() : null;
        String purpleDir = config.purpleDataDirectory();

        String qcFile = PurpleQCFile.generateFilename(purpleDir, tumorSample);
        String purityTsv = PurplePurity.generateFilename(purpleDir, tumorSample);
        String somaticDriverCatalogTsv = DriverCatalogFile.generateSomaticFilename(purpleDir, tumorSample);
        String somaticVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticVcfFile(purpleDir, tumorSample));
        String germlineDriverCatalogTsv = DriverCatalogFile.generateGermlineFilename(purpleDir, tumorSample);
        String germlineVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineVcfFile(purpleDir, tumorSample));
        String copyNumberTsv = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilename(purpleDir, tumorSample);
        String germlineDeletionTsv = GermlineDeletion.generateFilename(purpleDir, tumorSample);

        return load(
                tumorSample, referenceSample, rnaSample, qcFile, purityTsv,
                somaticDriverCatalogTsv, somaticVariantVcf, germlineDriverCatalogTsv, germlineVariantVcf,
                copyNumberTsv, geneCopyNumberTsv, germlineDeletionTsv, driverGenes);
    }

    private static String resolveVcfPath(final String vcfPath)
    {
        if(!new File(vcfPath).exists() && vcfPath.endsWith(".gz"))
        {
            String unzippedVcfPath = vcfPath.substring(0, vcfPath.length() - 3);
            if(new File(unzippedVcfPath).exists())
            {
                return unzippedVcfPath;
            }
        }
        return vcfPath;
    }

    private static PurpleData load(
            final String tumorSample, @Nullable String referenceSample, @Nullable String rnaSample,
            final String qcFile, final String purityTsv, final String somaticDriverCatalogTsv, final String somaticVariantVcf,
            final String germlineDriverCatalogTsv, final String germlineVariantVcf,
            final String copyNumberTsv, final String geneCopyNumberTsv,
            final String germlineDeletionTsv, final Map<String,DriverGene> driverGenes) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        List<PurpleVariantContext> allSomaticVariants = PurpleVariantContextLoader.withPassingOnlyFilter()
                .fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVariantVcf);

        List<PurpleVariantContext> panelSomaticVariants = allSomaticVariants.stream()
                .filter(x -> x.tier() == VariantTier.PANEL || x.tier() == VariantTier.HOTSPOT).collect(Collectors.toList());

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        geneCopyNumbers = geneCopyNumbers.stream().filter(x -> driverGenes.containsKey(x.GeneName)).collect(Collectors.toList());

        List<DriverCatalog> germlineDrivers = null;
        List<PurpleVariantContext> panelGermlineVariants = null;
        List<GermlineDeletion> panelGermlineDeletions = null;

        if(referenceSample != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);

            List<PurpleVariantContext> germlineVariants = new PurpleVariantContextLoader().fromVCFFile(
                    tumorSample, referenceSample, rnaSample, germlineVariantVcf);

            panelGermlineVariants = germlineVariants.stream()
                    .filter(x -> x.tier() == VariantTier.PANEL || x.tier() == VariantTier.HOTSPOT).collect(Collectors.toList());

            List<GermlineDeletion> germlineDeletions = GermlineDeletion.read(germlineDeletionTsv).stream()
                    .filter(x -> x.Filter.equals(CommonVcfTags.PASS_FILTER)).collect(Collectors.toList());

            panelGermlineDeletions = germlineDeletions.stream().filter(x -> driverGenes.containsKey(x.GeneName)).collect(Collectors.toList());
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .somaticVariants(panelSomaticVariants)
                .germlineVariants(panelGermlineVariants)
                .somaticCopyNumbers(PurpleCopyNumberFile.read(copyNumberTsv))
                .somaticGeneCopyNumbers(geneCopyNumbers)
                .germlineDeletions(panelGermlineDeletions)
                .build();
    }
}
