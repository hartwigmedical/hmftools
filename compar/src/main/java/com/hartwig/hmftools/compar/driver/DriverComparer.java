package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_PURPLE_SOMATIC_COPY_NUMBER;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.FieldCheckCache.getOrMakeFieldCheck;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_ABS_THRESHOLD;
import static com.hartwig.hmftools.compar.common.ComparConstants.COPY_NUMBER_PERC_THRESHOLD;
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.driver.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.common.purple.GeneCopyNumber;
import com.hartwig.hmftools.common.purple.GeneCopyNumberFile;
import com.hartwig.hmftools.common.purple.PurplePurity;
import com.hartwig.hmftools.common.purple.ReportedStatus;
import com.hartwig.hmftools.compar.common.CategoryType;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.PipelineSourcePaths;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

public class DriverComparer extends ItemComparer
{
    private final Map<SourceType,List<DriverCatalog>> mDrivers;
    private final Map<SourceType,PurplePurity> mPurities;
    private final Map<SourceType,List<GeneCopyNumber>> mGeneCopyNumbers;

    private boolean mPurpleSomaticDriversLoaded;
    private boolean mPurpleGermlineDriversLoaded;
    private boolean mLinxSomaticDriversLoaded;
    private boolean mLinxGermlineDriversLoaded;

    protected enum Fields
    {
        LikelihoodMethod,
        Likelihood,
        MinCopyNumber,
        MaxCopyNumber,
        Chromosome,
        ChromosomeBand,
        Purity,
        Ploidy;
    }

    public DriverComparer(final ComparConfig config, final Map<String, FieldCheck> fieldCheckMap)
    {
        super(config);

        mFields.add(new FieldInfo(
                Fields.LikelihoodMethod.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.LikelihoodMethod.toString()), null));

        mFields.add(new FieldInfo(
                Fields.Likelihood.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.Likelihood.toString(), 0.1, null),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MinCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MinCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.MaxCopyNumber.toString(),
                getOrMakeFieldCheck(fieldCheckMap, Fields.MaxCopyNumber.toString(), COPY_NUMBER_ABS_THRESHOLD, COPY_NUMBER_PERC_THRESHOLD),
                "%.2f"));

        mFields.add(new FieldInfo(
                Fields.Chromosome.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.Chromosome.toString()), null));

        mFields.add(new FieldInfo(
                Fields.ChromosomeBand.toString(), getOrMakeFieldCheck(fieldCheckMap, Fields.ChromosomeBand.toString()), null));

        mFields.add(FieldInfo.displayOnly(Fields.Purity.toString(), "%.2f"));
        mFields.add(FieldInfo.displayOnly(Fields.Ploidy.toString(), "%.2f"));

        mDrivers = Maps.newHashMap();
        mPurities = Maps.newHashMap();
        mGeneCopyNumbers = Maps.newHashMap();

        mPurpleSomaticDriversLoaded = false;
        mPurpleGermlineDriversLoaded = false;
        mLinxSomaticDriversLoaded = false;
        mLinxGermlineDriversLoaded = false;
    }

    @Override
    public CategoryType category() { return DRIVER; }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        // load data ahead of the standard calls for cross-driver comparisons
        for(SourceData sourceData : mConfig.Sources)
        {
            String sourceSampleId = mConfig.sourceSampleId(sourceData.Type, sampleId);

            String sourceReferenceId = mConfig.sourceReferenceId(sourceData.Type, sampleId);
            PipelineSourcePaths sampleFileSources = PipelineSourcePaths.sampleInstance(sourceData.PipelinePaths, sourceSampleId, sourceReferenceId);
            loadData(sourceSampleId, sampleFileSources, sourceData.Type);
        }

        boolean valid = CommonUtils.processSample(this, mConfig, sampleId, mismatches);

        mGeneCopyNumbers.clear();
        mPurities.clear();
        mDrivers.clear();
        return valid;
    }

    private void loadData(final String sampleId, final PipelineSourcePaths fileSources, final SourceType sourceType)
    {
        try
        {
            List<DriverCatalog> drivers = Lists.newArrayList();
            mDrivers.put(sourceType, drivers);

            String linxSomaticDriverFile = generateLinxSomaticDriverFilename(sampleId, fileSources);
            if(Files.exists(Paths.get(linxSomaticDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(linxSomaticDriverFile)
                        .stream()
                        .filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver()))
                        .collect(Collectors.toList()));
                mLinxSomaticDriversLoaded = true;
            }

            String purpleSomaticDriverFile = generatePurpleSomaticDriverFilename(sampleId, fileSources);
            if(Files.exists(Paths.get(purpleSomaticDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(purpleSomaticDriverFile));
                mPurpleSomaticDriversLoaded = true;
            }

            // add germline as well if present
            String purpleGermlineDriverFile = generatePurpleGermlineDriverFilename(sampleId, fileSources);
            if(Files.exists(Paths.get(purpleGermlineDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(purpleGermlineDriverFile));
                mPurpleGermlineDriversLoaded = true;
            }

            String linxGermlineDriverFile = generateLinxGermlineDriverFilename(sampleId, fileSources);
            if(Files.exists(Paths.get(linxGermlineDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(linxGermlineDriverFile)
                        .stream()
                        .filter(x -> DRIVERS_LINX_GERMLINE.contains(x.driver()))
                        .collect(Collectors.toList()));
                mLinxGermlineDriversLoaded = true;
            }

            PurplePurity purity = PurplePurity.read(PurplePurity.generateFilename(fileSources.Purple, sampleId));
            mPurities.put(sourceType, purity);

            String geneCopyNumberFile = GeneCopyNumberFile.generateFilename(fileSources.Purple, sampleId);
            mGeneCopyNumbers.put(sourceType, GeneCopyNumberFile.read(geneCopyNumberFile));
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load driver data: {}", sampleId, e.toString());
        }
    }

    @Override
    public List<String> displayFieldNames()
    {
        return Arrays.stream(Fields.values()).map(x -> x.toString()).collect(Collectors.toList());
    }

    private List<ComparableItem> createDriverItems(final SourceType sourceType)
    {
        List<DriverCatalog> drivers = mDrivers.get(sourceType);

        if(drivers == null)
            return null;

        List<ComparableItem> items = Lists.newArrayList();

        for(DriverCatalog driverCatalog : drivers)
        {
            items.add(createDriverData(driverCatalog, mPurities.get(sourceType), true));
        }

        // create non-reportable CN driver events if present in the other source
        List<DriverCatalog> otherDrivers = mDrivers.get(sourceType == OLD ? NEW : OLD);

        for(DriverCatalog otherDriver : otherDrivers)
        {
            if(otherDriver.reportedStatus() != ReportedStatus.REPORTED)
                continue;

            if(!DRIVERS_PURPLE_SOMATIC_COPY_NUMBER.contains(otherDriver.driver()))
                continue;

            if(drivers.stream().anyMatch(x -> x.driver() == otherDriver.driver() && x.gene().equals(otherDriver.gene())))
                continue;

            GeneCopyNumber geneCopyNumber = mGeneCopyNumbers.get(sourceType).stream()
                    .filter(x -> x.GeneName.equals(otherDriver.gene())).findFirst().orElse(null);

            if(geneCopyNumber == null)
                continue;

            ImmutableDriverCatalog.Builder builder = ImmutableDriverCatalog.builder().from(otherDriver);
            builder.minCopyNumber(geneCopyNumber.MinCopyNumber);
            builder.maxCopyNumber(geneCopyNumber.MaxCopyNumber);
            builder.driverLikelihood(0);
            builder.reportedStatus(ReportedStatus.NOT_REPORTED);

            items.add(createDriverData(builder.build(), mPurities.get(sourceType), false));
        }

        return items;
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final PipelineSourcePaths fileSources)
    {
        Set<CategoryType> categories = mConfig.Categories;

        boolean missingCatalogFile = false;
        if(!mPurpleSomaticDriversLoaded && hasOverlap(categories, CategoryType.purpleSomaticOnlyCategories()))
        {
            CMP_LOGGER.warn(
                    "sample({}) failed to find and load Purple somatic driver catalog file: {}", sampleId,
                    generatePurpleSomaticDriverFilename(sampleId, fileSources));
            missingCatalogFile = true;
        }

        if(!mLinxSomaticDriversLoaded && hasOverlap(categories, CategoryType.linxSomaticOnlyCategories()))
        {
            CMP_LOGGER.warn(
                    "sample({}) failed to find and load Linx somatic driver catalog file: {}", sampleId,
                    generateLinxSomaticDriverFilename(sampleId, fileSources));
            missingCatalogFile = true;
        }

        if(!mPurpleGermlineDriversLoaded && hasOverlap(categories, CategoryType.purpleGermlineOnlyCategories()))
        {
            CMP_LOGGER.warn(
                    "sample({}) failed to find and load Purple germline driver catalog file: {}", sampleId,
                    generatePurpleGermlineDriverFilename(sampleId, fileSources));
            missingCatalogFile = true;
        }

        if(!mLinxGermlineDriversLoaded && hasOverlap(categories, CategoryType.linxGermlineOnlyCategories()))
        {
            CMP_LOGGER.warn(
                    "sample({}) failed to find and load Linx germline driver catalog file: {}", sampleId,
                    generateLinxGermlineDriverFilename(sampleId, fileSources));
            missingCatalogFile = true;
        }

        if(missingCatalogFile)
        {
            return null;
        }
        else
        {
            return createDriverItems(fileSources.Source);
        }
    }

    private DriverData createDriverData(final DriverCatalog driver, final PurplePurity purity, final boolean isPass)
    {
        boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
        String comparisonChromosome = determineComparisonChromosome(driver.chromosome(), mConfig.RequiresLiftover);
        return new DriverData(driver, purity, comparisonChromosome, checkTranscript, isPass, mFields);
    }

    private static String generatePurpleSomaticDriverFilename(final String sampleId, final PipelineSourcePaths fileSources)
    {
        return DriverCatalogFile.generateSomaticFilename(fileSources.Purple, sampleId);
    }

    private static String generateLinxSomaticDriverFilename(final String sampleId, final PipelineSourcePaths fileSources)
    {
        return LinxDriver.generateCatalogFilename(fileSources.Linx, sampleId, true);
    }

    private static String generatePurpleGermlineDriverFilename(final String sampleId, final PipelineSourcePaths fileSources)
    {
        return DriverCatalogFile.generateGermlineFilename(fileSources.Purple, sampleId);
    }

    private static String generateLinxGermlineDriverFilename(final String sampleId, final PipelineSourcePaths fileSources)
    {
        return LinxDriver.generateCatalogFilename(fileSources.LinxGermline, sampleId, false);
    }

    private static boolean hasOverlap(final Set<CategoryType> set1, final Set<CategoryType> set2)
    {
        return !Sets.intersection(set1, set2).isEmpty();
    }
}
