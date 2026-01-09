package com.hartwig.hmftools.orange.algo.purple;

import static com.hartwig.hmftools.common.sv.SvVcfTags.INFERRED;

import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.panel.DriverGene;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.GermlineAmpDel;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.purple.PurpleCommon;
import com.hartwig.hmftools.common.purple.PurpleCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.PurpleQCFile;
import com.hartwig.hmftools.common.sv.StructuralVariant;
import com.hartwig.hmftools.common.sv.StructuralVariantFileLoader;
import com.hartwig.hmftools.common.variant.CommonVcfTags;
import com.hartwig.hmftools.common.variant.VariantTier;
import com.hartwig.hmftools.orange.OrangeConfig;

import org.jetbrains.annotations.Nullable;

import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.variantcontext.filter.PassingVariantFilter;

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
        String somaticStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleSomaticSvFile(purpleDir, tumorSample));
        String germlineStructuralVariantVcf = resolveVcfPath(PurpleCommon.purpleGermlineSvFile(purpleDir, tumorSample));
        String copyNumberTsv = PurpleCopyNumberFile.generateFilenameForReading(purpleDir, tumorSample);
        String geneCopyNumberTsv = GeneCopyNumberFile.generateFilename(purpleDir, tumorSample);
        String germlineDeletionTsv = GermlineAmpDel.generateFilename(purpleDir, tumorSample);
        String segmentTsv = SegmentFile.generateFilename(purpleDir, tumorSample);

        return load(
                tumorSample, referenceSample, rnaSample,
                qcFile, purityTsv, somaticDriverCatalogTsv, somaticVariantVcf, germlineDriverCatalogTsv, germlineVariantVcf,
                somaticStructuralVariantVcf, germlineStructuralVariantVcf,
                copyNumberTsv, geneCopyNumberTsv, germlineDeletionTsv, segmentTsv,
                driverGenes, config.includeNonGenePanelEvents());
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
            final String germlineDriverCatalogTsv, final String germlineVariantVcf, final String somaticStructuralVariantVcf,
            final String germlineStructuralVariantVcf, final String copyNumberTsv, final String geneCopyNumberTsv,
            final String germlineDeletionTsv, final String segmentTsv,
            final Map<String,DriverGene> driverGenes, final boolean includeNonGenePanelEvents) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(qcFile, purityTsv);

        List<DriverCatalog> somaticDrivers = DriverCatalogFile.read(somaticDriverCatalogTsv);

        List<PurpleVariantContext> allSomaticVariants = PurpleVariantContextLoader.withPassingOnlyFilter()
                .fromVCFFile(tumorSample, referenceSample, rnaSample, somaticVariantVcf);

        List<PurpleVariantContext> panelSomaticVariants = allSomaticVariants.stream()
                .filter(x -> x.tier() == VariantTier.PANEL || x.tier() == VariantTier.HOTSPOT).collect(Collectors.toList());

        if(!includeNonGenePanelEvents)
            allSomaticVariants.clear();

        List<StructuralVariant> somaticPassingSVs, somaticInferredSVs;

        if(includeNonGenePanelEvents)
        {
            somaticPassingSVs = Lists.newArrayList();
            somaticInferredSVs = Lists.newArrayList();
            loadStructuralVariants(somaticStructuralVariantVcf, somaticPassingSVs, somaticInferredSVs);
        }
        else
        {
            somaticPassingSVs = Collections.emptyList();
            somaticInferredSVs = Collections.emptyList();
        }

        List<GeneCopyNumber> geneCopyNumbers = GeneCopyNumberFile.read(geneCopyNumberTsv);

        if(!includeNonGenePanelEvents)
            geneCopyNumbers = geneCopyNumbers.stream().filter(x -> driverGenes.containsKey(x.GeneName)).collect(Collectors.toList());

        List<Segment> segments = includeNonGenePanelEvents ? SegmentFile.read(segmentTsv) : Collections.emptyList();

        List<DriverCatalog> germlineDrivers = null;
        List<PurpleVariantContext> allGermlineVariants = null;
        List<PurpleVariantContext> panelGermlineVariants = null;
        List<GermlineAmpDel> allGermlineDeletions = null;
        List<GermlineAmpDel> panelGermlineDeletions = null;

        List<StructuralVariant> germlinePassingSVs = null;
        List<StructuralVariant> germlineInferredSVs = null;

        if(referenceSample != null)
        {
            germlineDrivers = DriverCatalogFile.read(germlineDriverCatalogTsv);

            allGermlineVariants = new PurpleVariantContextLoader().fromVCFFile(tumorSample, referenceSample, rnaSample, germlineVariantVcf);

            panelGermlineVariants = allGermlineVariants.stream()
                    .filter(x -> x.tier() == VariantTier.PANEL || x.tier() == VariantTier.HOTSPOT).collect(Collectors.toList());

            if(!includeNonGenePanelEvents)
                allGermlineVariants.clear();

            allGermlineDeletions = GermlineAmpDel.read(germlineDeletionTsv).stream()
                    .filter(x -> x.Filter.equals(CommonVcfTags.PASS_FILTER)).collect(Collectors.toList());

            panelGermlineDeletions = allGermlineDeletions.stream().filter(x -> driverGenes.containsKey(x.GeneName)).collect(Collectors.toList());

            if(!includeNonGenePanelEvents)
                allGermlineDeletions.clear();

            if(includeNonGenePanelEvents)
            {
                germlinePassingSVs = Lists.newArrayList();
                germlineInferredSVs = Lists.newArrayList();
                loadStructuralVariants(germlineStructuralVariantVcf, somaticPassingSVs, somaticInferredSVs);
            }
            else
            {
                germlinePassingSVs = Collections.emptyList();
                germlineInferredSVs = Collections.emptyList();
            }
        }

        return ImmutablePurpleData.builder()
                .purityContext(purityContext)
                .somaticDrivers(somaticDrivers)
                .germlineDrivers(germlineDrivers)
                .allSomaticVariants(allSomaticVariants)
                .driverSomaticVariants(panelSomaticVariants)
                .allGermlineVariants(allGermlineVariants)
                .driverGermlineVariants(panelGermlineVariants)
                .allPassingSomaticStructuralVariants(somaticPassingSVs)
                .allPassingGermlineStructuralVariants(germlinePassingSVs)
                .allInferredSomaticStructuralVariants(somaticInferredSVs)
                .allInferredGermlineStructuralVariants(germlineInferredSVs)
                .somaticCopyNumbers(PurpleCopyNumberFile.read(copyNumberTsv))
                .somaticGeneCopyNumbers(geneCopyNumbers)
                .allGermlineDeletions(allGermlineDeletions)
                .driverGermlineDeletions(panelGermlineDeletions)
                .segments(segments)
                .build();
    }

    private static void loadStructuralVariants(
            final String vcfPath, final List<StructuralVariant> passingSVs, final List<StructuralVariant> inferredSVs) throws IOException
    {
        CompoundFilter passingOrInferredFilter = new CompoundFilter(false);
        passingOrInferredFilter.add(new PassingVariantFilter());
        passingOrInferredFilter.add(variantContext -> variantContext.getFilters().contains(INFERRED));

        List<StructuralVariant> passingOrInferred = StructuralVariantFileLoader.fromFile(vcfPath, passingOrInferredFilter);

        for(StructuralVariant variant : passingOrInferred)
        {
            if(variant.filter().equals(INFERRED))
                inferredSVs.add(variant);
            else
                passingSVs.add(variant);
        }
    }
}
