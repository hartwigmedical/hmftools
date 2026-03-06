package com.hartwig.hmftools.orange.algo.purple;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.ChrArmCopyNumber;
import com.hartwig.hmftools.common.purple.ChrArmCopyNumbersFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.common.variant.SmallVariant;
import com.hartwig.hmftools.common.variant.SmallVariantFactory;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.Nullable;

public final class PurpleDataLoader
{
    public static PurpleData load(final OrangeConfig config, final Map<String,DriverGene> driverGenes) throws IOException
    {
        String tumorId = config.TumorId;
        String purpleDir = config.PurpleDataDirectory;

        String qcFile = PurpleQCFile.generateFilename(purpleDir, tumorId);
        String purityTsv = PurplePurity.generateFilename(purpleDir, tumorId);
        String somaticDriverCatalogTsv = DriverCatalogFile.generateSomaticFilename(purpleDir, tumorId);
        String somaticVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticVcfFile(purpleDir, tumorId));
        String germlineDriverCatalogTsv = DriverCatalogFile.generateGermlineFilename(purpleDir, tumorId);
        String germlineVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineVcfFile(purpleDir, tumorId));
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilename(purpleDir, tumorId);
        String germlineDeletionTsv = GermlineAmpDel.generateFilename(purpleDir, tumorId);

        String chrArmCopyNumberTsv = ChrArmCopyNumbersFile.generateFilename(purpleDir, tumorId);

        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        // exclude non-reportable events
        somaticDrivers = somaticDrivers.stream().filter(x -> x.reportedStatus() != ReportedStatus.NOT_REPORTED).collect(Collectors.toList());

        List<SmallVariant> allSomaticVariants = SmallVariantFactory.passOnlyInstance().fromVCFFile(
                tumorId, config.ReferenceId, config.RnaSampleId, somaticVariantVcf);

        List<SmallVariant> panelSomaticVariants = allSomaticVariants.stream().filter(x -> x.reported()).collect(Collectors.toList());

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        geneCopyNumbers = geneCopyNumbers.stream().filter(x -> driverGenes.containsKey(x.GeneName)).collect(Collectors.toList());

        List<ChrArmCopyNumber> chrArmCopyNumbers = ChrArmCopyNumbersFile.read(chrArmCopyNumberTsv);

        List<DriverCatalog> germlineDrivers = null;
        List<SmallVariant> panelGermlineVariants = null;
        List<GermlineAmpDel> panelGermlineAmpDels = null;

        if(config.hasReference())
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);

            germlineDrivers = germlineDrivers.stream().filter(x -> x.reportedStatus() != ReportedStatus.NOT_REPORTED).collect(Collectors.toList());

            List<SmallVariant> germlineVariants = new SmallVariantFactory().fromVCFFile(
                    tumorId, config.ReferenceId, config.RnaSampleId, germlineVariantVcf);

            panelGermlineVariants = germlineVariants.stream().filter(x -> x.reported()).collect(Collectors.toList());

            List<GermlineAmpDel> germlineAmpDels = GermlineAmpDel.read(germlineDeletionTsv).stream()
                    .filter(x -> x.Filter.equals(CommonVcfTags.PASS_FILTER)).collect(Collectors.toList());

            panelGermlineAmpDels = germlineAmpDels.stream().filter(x -> x.Reported == ReportedStatus.REPORTED).collect(Collectors.toList());
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .somaticVariants(panelSomaticVariants)
                .germlineVariants(panelGermlineVariants)
                .somaticGeneCopyNumbers(geneCopyNumbers)
                .germlineAmpDels(panelGermlineAmpDels)
                .chrArmCopyNumbers(chrArmCopyNumbers)
                .build();
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
}
