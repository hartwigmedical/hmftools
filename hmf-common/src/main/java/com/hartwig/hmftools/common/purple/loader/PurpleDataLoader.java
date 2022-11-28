package com.hartwig.hmftools.common.purple.loader;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

public final class PurpleDataLoader
{
    private PurpleDataLoader() {}

    @NotNull
    public static PurpleData load(final String tumorSample, @Nullable final String referenceSample, @Nullable final String rnaSample,
            final String purpleDir) throws IOException
    {
        String qcFile = PurpleQCFile.generateFilename(purpleDir, tumorSample);
        String purityTsv = PurityContextFile.generateFilenameForReading(purpleDir, tumorSample);
        String somaticDriverCatalogTsv = DriverCatalogFile.generateSomaticFilename(purpleDir, tumorSample);
        String somaticVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticVcfFile(purpleDir, tumorSample));
        String germlineDriverCatalogTsv = DriverCatalogFile.generateGermlineFilename(purpleDir, tumorSample);
        String germlineVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineVcfFile(purpleDir, tumorSample));
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String copyNumberTsv = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String germlineDeletionTsv = GermlineDeletion.generateFilename(purpleDir, tumorSample);

        return load(tumorSample, referenceSample, rnaSample, qcFile, purityTsv, somaticDriverCatalogTsv, somaticVariantVcf,
                germlineDriverCatalogTsv, germlineVariantVcf, geneCopyNumberTsv, copyNumberTsv, germlineDeletionTsv);
    }

    private static String resolveVcfPath(final String vcfPath)
    {
        if (!new File(vcfPath).exists() && vcfPath.endsWith(".gz"))
        {
            String unzippedVcfPath = vcfPath.substring(0, vcfPath.length() - 3);
            if (new File(unzippedVcfPath).exists())
            {
                return unzippedVcfPath;
            }
        }
        return vcfPath;
    }

    @NotNull
    private static PurpleData load(@NotNull String tumorSample, @Nullable String referenceSample, @Nullable String rnaSample,
            @NotNull String qcFile, @NotNull String purityTsv, @NotNull String somaticDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String germlineDriverCatalogTsv, @NotNull String germlineVariantVcf, @NotNull String geneCopyNumberTsv,
            @NotNull String copyNumberTsv, @NotNull String germlineDeletionTsv) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        List<SomaticVariant> allSomaticVariants =
                SomaticVariantFactory.passOnlyInstance().fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVariantVcf);
        List<SomaticVariant> reportableSomaticVariants = selectReportedVariants(allSomaticVariants);

        List<GeneCopyNumber> allSomaticGeneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        List<DriverCatalog> germlineDrivers = null;
        List<SomaticVariant> allGermlineVariants = null;
        List<SomaticVariant> reportableGermlineVariants = null;
        List<GermlineDeletion> allGermlineDeletions = null;
        List<GermlineDeletion> reportableGermlineDeletions = null;
        if(referenceSample != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);

            allGermlineVariants = new SomaticVariantFactory().fromVCFFile(tumorSample, referenceSample, rnaSample, germlineVariantVcf);
            reportableGermlineVariants = selectReportedVariants(allGermlineVariants);

            allGermlineDeletions = GermlineDeletion.read(germlineDeletionTsv);
            reportableGermlineDeletions = selectReportedDeletions(allGermlineDeletions);
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .allSomaticVariants(allSomaticVariants)
                .reportableSomaticVariants(reportableSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .allSomaticCopyNumbers(PurpleCopyNumberFile.read(copyNumberTsv))
                .allSomaticGeneCopyNumbers(allSomaticGeneCopyNumbers)
                .allGermlineDeletions(allGermlineDeletions)
                .reportableGermlineDeletions(reportableGermlineDeletions)
                .build();
    }

    @NotNull
    private static List<SomaticVariant> selectReportedVariants(@NotNull List<SomaticVariant> allVariants)
    {
        List<SomaticVariant> reported = Lists.newArrayList();
        for(SomaticVariant variant : allVariants)
        {
            if(variant.reported())
            {
                reported.add(variant);
            }
        }
        return reported;
    }

    @NotNull
    private static List<GermlineDeletion> selectReportedDeletions(@NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        List<GermlineDeletion> reported = Lists.newArrayList();
        for(GermlineDeletion deletion : allGermlineDeletions)
        {
            if(deletion.Reported)
            {
                reported.add(deletion);
            }
        }
        return reported;
    }
}
