package com.hartwig.hmftools.common.purple.loader;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.drivercatalog.AmplificationDrivers;
import com.hartwig.hmftools.common.drivercatalog.DeletionDrivers;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogKey;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogMap;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGeneGermlineReporting;
import com.hartwig.hmftools.common.drivercatalog.panel.DriverGenePanel;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGene;
import com.hartwig.hmftools.common.drivercatalog.panel.ImmutableDriverGenePanel;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.purple.PurpleQCStatus;
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
        List<GainLoss> allSomaticGainsLosses = extractAllGainsLosses(
                    purityContext.qc().status(), purityContext.bestFit().ploidy(), purityContext.targeted(), allSomaticGeneCopyNumbers);
        List<GainLoss> reportableSomaticGainsLosses = somaticGainsLossesFromDrivers(somaticDrivers);

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
                .allSomaticGainsLosses(allSomaticGainsLosses)
                .reportableSomaticGainsLosses(reportableSomaticGainsLosses)
                .allGermlineDeletions(allGermlineDeletions)
                .reportableGermlineDeletions(reportableGermlineDeletions)
                .build();
    }

    @NotNull
    private static List<GainLoss> somaticGainsLossesFromDrivers(@NotNull List<DriverCatalog> drivers)
    {
        List<GainLoss> gainsLosses = Lists.newArrayList();

        Map<DriverCatalogKey, DriverCatalog> geneDriverMap = DriverCatalogMap.toDriverMap(drivers);
        for(DriverCatalogKey key : geneDriverMap.keySet())
        {
            DriverCatalog geneDriver = geneDriverMap.get(key);

            if(geneDriver.driver() == DriverType.AMP || geneDriver.driver() == DriverType.PARTIAL_AMP
                    || geneDriver.driver() == DriverType.DEL)
            {
                gainsLosses.add(toGainLoss(geneDriver));
            }
        }
        return gainsLosses;
    }

    private static List<GainLoss> extractAllGainsLosses(final Set<PurpleQCStatus> qcStatus, double ploidy, boolean isTargetRegions,
            final List<GeneCopyNumber> geneCopyNumbers)
    {
        List<DriverGene> allGenes = Lists.newArrayList();
        for(GeneCopyNumber geneCopyNumber : geneCopyNumbers)
        {
            allGenes.add(ImmutableDriverGene.builder()
                    .gene(geneCopyNumber.geneName())
                    .reportMissenseAndInframe(false)
                    .reportNonsenseAndFrameshift(false)
                    .reportSplice(false)
                    .reportDeletion(true)
                    .reportDisruption(false)
                    .reportAmplification(true)
                    .reportSomaticHotspot(false)
                    .reportGermlineVariant(DriverGeneGermlineReporting.NONE)
                    .reportGermlineHotspot(DriverGeneGermlineReporting.NONE)
                    .reportGermlineDisruption(false)
                    .likelihoodType(DriverCategory.ONCO)
                    .reportPGX(false)
                    .build());
        }

        final DriverGenePanel allGenesPanel = ImmutableDriverGenePanel.builder().driverGenes(allGenes).build();
        AmplificationDrivers ampDrivers = new AmplificationDrivers(qcStatus, allGenesPanel);
        DeletionDrivers delDrivers = new DeletionDrivers(qcStatus, allGenesPanel);

        List<DriverCatalog> allGainLosses = Lists.newArrayList();
        allGainLosses.addAll(ampDrivers.amplifications(ploidy, geneCopyNumbers, isTargetRegions));
        allGainLosses.addAll(delDrivers.deletions(geneCopyNumbers, isTargetRegions));

        return somaticGainsLossesFromDrivers(allGainLosses);
    }

    @NotNull
    private static GainLoss toGainLoss(@NotNull DriverCatalog driver)
    {
        return ImmutableGainLoss.builder()
                .chromosome(driver.chromosome())
                .chromosomeBand(driver.chromosomeBand())
                .gene(driver.gene())
                .transcript(driver.transcript())
                .isCanonical(driver.isCanonical())
                .interpretation(CopyNumberInterpretation.fromCNADriver(driver))
                .minCopies(Math.round(Math.max(0, driver.minCopyNumber())))
                .maxCopies(Math.round(Math.max(0, driver.maxCopyNumber())))
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
