package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Objects;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineDeletion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;

import org.jetbrains.annotations.NotNull;
import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

public final class PurpleDataLoader
{
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
        String somaticStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticSvFile(purpleDir, tumorSample));
        String germlineStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineSvFile(purpleDir, tumorSample));
        String copyNumberTsv = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilename(purpleDir, tumorSample);
        String germlineDeletionTsv = GermlineDeletion.generateFilename(purpleDir, tumorSample);
        String segmentTsv = SegmentFile.generateFilename(purpleDir, tumorSample);

        return load(tumorSample,
                referenceSample,
                rnaSample,
                qcFile,
                purityTsv,
                somaticDriverCatalogTsv,
                somaticVariantVcf,
                germlineDriverCatalogTsv,
                germlineVariantVcf,
                somaticStructuralVariantVcf,
                germlineStructuralVariantVcf,
                copyNumberTsv,
                geneCopyNumberTsv,
                germlineDeletionTsv,
                segmentTsv);
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

    @NotNull
    private static PurpleData load(@NotNull String tumorSample, @Nullable String referenceSample, @Nullable String rnaSample,
            @NotNull String qcFile, @NotNull String purityTsv, @NotNull String somaticDriverCatalogTsv, @NotNull String somaticVariantVcf,
            @NotNull String germlineDriverCatalogTsv, @NotNull String germlineVariantVcf, @NotNull String somaticStructuralVariantVcf,
            @NotNull String germlineStructuralVariantVcf, @NotNull String copyNumberTsv, @NotNull String geneCopyNumberTsv,
            @NotNull String germlineDeletionTsv, @NotNull String segmentTsv) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        List<PurpleVariantContext> allSomaticVariants = PurpleVariantContextLoader.withPassingOnlyFilter()
                .fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVariantVcf);
        List<PurpleVariantContext> reportableSomaticVariants = selectReportedVariants(allSomaticVariants);
        StructuralVariants somaticStructuralVariants = loadStructuralVariants(somaticStructuralVariantVcf);

        List<GeneCopyNumber> allSomaticGeneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        List<Segment> segments = SegmentFile.read(segmentTsv);

        List<DriverCatalog> germlineDrivers = null;
        List<PurpleVariantContext> allGermlineVariants = null;
        List<PurpleVariantContext> reportableGermlineVariants = null;
        List<GermlineDeletion> allGermlineDeletions = null;
        List<GermlineDeletion> reportableGermlineDeletions = null;
        StructuralVariants germlineStructuralVariants = new StructuralVariants(null, null);
        if(referenceSample != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);

            allGermlineVariants = new PurpleVariantContextLoader().fromVCFFile(tumorSample, referenceSample, rnaSample,
                    germlineVariantVcf);
            reportableGermlineVariants = selectReportedVariants(allGermlineVariants);

            allGermlineDeletions = selectPassDeletions(GermlineDeletion.read(germlineDeletionTsv));
            reportableGermlineDeletions = selectReportedDeletions(allGermlineDeletions);

            germlineStructuralVariants = loadStructuralVariants(germlineStructuralVariantVcf);
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .allSomaticVariants(allSomaticVariants)
                .reportableSomaticVariants(reportableSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .reportableGermlineVariants(reportableGermlineVariants)
                .allPassingSomaticStructuralVariants(somaticStructuralVariants.allPassingStructuralVariants)
                .allPassingGermlineStructuralVariants(germlineStructuralVariants.allPassingStructuralVariants)
                .allInferredSomaticStructuralVariants(somaticStructuralVariants.allInferredStructuralVariants)
                .allInferredGermlineStructuralVariants(germlineStructuralVariants.allInferredStructuralVariants)
                .allSomaticCopyNumbers(PurpleCopyNumberFile.read(copyNumberTsv))
                .allSomaticGeneCopyNumbers(allSomaticGeneCopyNumbers)
                .allGermlineDeletions(allGermlineDeletions)
                .reportableGermlineDeletions(reportableGermlineDeletions)
                .segments(segments)
                .build();
    }

    @NotNull
    private static List<PurpleVariantContext> selectReportedVariants(@NotNull List<PurpleVariantContext> allVariants)
    {
        List<PurpleVariantContext> reported = Lists.newArrayList();
        for(PurpleVariantContext variant : allVariants)
        {
            if(variant.reported())
            {
                reported.add(variant);
            }
        }
        return reported;
    }

    @NotNull
    private static List<GermlineDeletion> selectPassDeletions(@NotNull List<GermlineDeletion> allGermlineDeletions)
    {
        List<GermlineDeletion> pass = Lists.newArrayList();
        for(GermlineDeletion deletion : allGermlineDeletions)
        {
            if(deletion.Filter.equals("PASS"))
            {
                pass.add(deletion);
            }
        }

        return pass;
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

    private static class StructuralVariants
    {
        @Nullable
        final List<StructuralVariant> allPassingStructuralVariants;
        @Nullable
        final List<StructuralVariant> allInferredStructuralVariants;

        private StructuralVariants(@Nullable final List<StructuralVariant> allPassingStructuralVariants,
                @Nullable List<StructuralVariant> allInferredStructuralVariants)
        {
            this.allPassingStructuralVariants = allPassingStructuralVariants;
            this.allInferredStructuralVariants = allInferredStructuralVariants;
        }
    }

    @NotNull
    private static StructuralVariants loadStructuralVariants(@NotNull String vcfPath) throws IOException
    {
        CompoundFilter passingOrInferredFilter = new CompoundFilter(false);
        passingOrInferredFilter.add(new PassingVariantFilter());
        passingOrInferredFilter.add(variantContext -> variantContext.getFilters().contains(INFERRED));

        List<StructuralVariant> passingOrInferred =
                StructuralVariantFileLoader.fromFile(vcfPath, passingOrInferredFilter);

        Map<Boolean, List<StructuralVariant>> partitioned =
                passingOrInferred.stream().collect(Collectors.partitioningBy(x -> Objects.equals(x.filter(), INFERRED)));

        return new StructuralVariants(partitioned.get(false), partitioned.get(true));
    }
}
