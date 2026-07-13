package com.hartwig.hmftools.compar.driver;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_PURPLE_SOMATIC_COPY_NUMBER;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.compar.common.CategoryType.DRIVER;
import static com.hartwig.hmftools.compar.common.CommonUtils.FLD_CHROMOSOME_BAND;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.common.SourceType.NEW;
import static com.hartwig.hmftools.compar.common.SourceType.OLD;
import static com.hartwig.hmftools.compar.purple.PurityComparer.FLD_PLOIDY;
import static com.hartwig.hmftools.compar.purple.PurityComparer.FLD_PURITY;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
import com.hartwig.hmftools.compar.common.FieldConfig;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.compar.common.SourceData;
import com.hartwig.hmftools.compar.common.SourceType;
import com.hartwig.hmftools.compar.common.field.DisplayField;
import com.hartwig.hmftools.compar.common.field.DoubleField;
import com.hartwig.hmftools.compar.common.field.Field;
import com.hartwig.hmftools.compar.common.field.StringField;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverComparer implements ItemComparer
{
    protected static final String FLD_LIKELIHOOD = "Likelihood";
    protected static final String FLD_LIKE_METHOD = "LikelihoodMethod";
    protected static final String FLD_MIN_COPY_NUMBER = "MinCopyNumber";
    protected static final String FLD_MAX_COPY_NUMBER = "MaxCopyNumber";

    private final ComparConfig mConfig;

    private final Map<SourceType,List<DriverCatalog>> mDrivers;
    private final Map<SourceType,PurplePurity> mPurities;
    private final Map<SourceType,List<GeneCopyNumber>> mGeneCopyNumbers;

    public DriverComparer(final ComparConfig config)
    {
        mConfig = config;

        mDrivers = Maps.newHashMap();
        mPurities = Maps.newHashMap();
        mGeneCopyNumbers = Maps.newHashMap();
    }

    @Override
    public CategoryType category() { return DRIVER; }

    @Override
    public List<Field> fields(final MatchLevel matchLevel)
    {
        return List.of(
                new StringField(FLD_LIKE_METHOD, i -> ((DriverData) i).DriverCatalog.likelihoodMethod().toString(),
                        true),
                new DoubleField(FLD_LIKELIHOOD, i -> ((DriverData) i).DriverCatalog.driverLikelihood(), true,
                        0.1, null, "%.2f"),
                new DoubleField(FLD_MIN_COPY_NUMBER, i -> ((DriverData) i).DriverCatalog.minCopyNumber(), true,
                        0.3, 0.15, "%.2f"),
                new DoubleField(FLD_MAX_COPY_NUMBER, i -> ((DriverData) i).DriverCatalog.maxCopyNumber(), true,
                        0.3, 0.15, "%.2f"),
                new StringField(FLD_CHROMOSOME, i -> ((DriverData) i).mComparisonChromosome, true),
                new StringField(FLD_CHROMOSOME_BAND, i -> ((DriverData) i).DriverCatalog.chromosomeBand(), true),
                new DisplayField(FLD_PURITY, i -> format("%.2f", ((DriverData) i).mPurity.Purity), i -> true),
                new DisplayField(FLD_PLOIDY, i -> format("%.2f", ((DriverData) i).mPurity.Ploidy), i -> true)
        );
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches, final FieldConfig fieldConfig)
    {
        // load data ahead of the standard calls for cross-driver comparisons
        for(SourceData sourceData : mConfig.Sources)
        {
            String sourceSampleId = mConfig.sourceSampleId(sourceData.Type, sampleId);

            if(sourceData.Database != null)
            {
                loadData(sourceSampleId, sourceData.Database, sourceData.Type);
            }
            else
            {
                String sourceReferenceId = mConfig.sourceReferenceId(sourceData.Type, sampleId);
                FileSources sampleFileSources = FileSources.sampleInstance(sourceData.Files, sourceSampleId, sourceReferenceId);
                loadData(sourceSampleId, sampleFileSources, sourceData.Type);
            }
        }

        boolean valid = CommonUtils.processSample(this, mConfig, sampleId, mismatches, fieldConfig);

        mGeneCopyNumbers.clear();
        mPurities.clear();
        mDrivers.clear();
        return valid;
    }

    private void loadData(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        mDrivers.put(sourceType, Lists.newArrayList(dbAccess.readDriverCatalog(sampleId)));
    }

    private void loadData(final String sampleId, final FileSources fileSources, final SourceType sourceType)
    {
        try
        {
            // use Linx if present, otherwise Purple drivers
            String linxDriverFile = LinxDriver.generateCatalogFilename(fileSources.Linx, sampleId, true);
            String purpleDriverFile = DriverCatalogFile.generateSomaticFilename(fileSources.Purple, sampleId);

            List<DriverCatalog> drivers = Lists.newArrayList();
            mDrivers.put(sourceType, drivers);

            if(Files.exists(Paths.get(linxDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(linxDriverFile)
                        .stream()
                        .filter(x -> DRIVERS_LINX_SOMATIC.contains(x.driver()))
                        .collect(Collectors.toList()));
            }

            if(Files.exists(Paths.get(purpleDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(purpleDriverFile));
            }

            // add germline as well if present
            String purpleGermlineDriverFile = DriverCatalogFile.generateGermlineFilename(fileSources.Purple, sampleId);
            String linxGermlineDriverFile = LinxDriver.generateCatalogFilename(fileSources.LinxGermline, sampleId, false);

            if(Files.exists(Paths.get(purpleGermlineDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(purpleGermlineDriverFile));
            }

            if(Files.exists(Paths.get(linxGermlineDriverFile)))
            {
                drivers.addAll(DriverCatalogFile.read(linxGermlineDriverFile)
                        .stream()
                        .filter(x -> DRIVERS_LINX_GERMLINE.contains(x.driver()))
                        .collect(Collectors.toList()));
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
        return Lists.newArrayList(FLD_LIKE_METHOD, FLD_LIKELIHOOD, FLD_MIN_COPY_NUMBER, FLD_MAX_COPY_NUMBER);
    }

    private List<ComparableItem> createDriverItems(final SourceType sourceType)
    {
        List<DriverCatalog> drivers = mDrivers.get(sourceType);

        if(drivers == null)
            return null;

        List<ComparableItem> items = Lists.newArrayList();

        for(DriverCatalog driverCatalog : drivers)
        {
            items.add(createDriverData(driverCatalog, mPurities.get(sourceType)));
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

            items.add(createDriverData(builder.build(), mPurities.get(sourceType)));
        }

        return items;
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final SourceType sourceType)
    {
        return createDriverItems(sourceType);
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        return createDriverItems(fileSources.Source);
    }

    private DriverData createDriverData(final DriverCatalog driver, final PurplePurity purity)
    {
        boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
        String comparisonChromosome = determineComparisonChromosome(driver.chromosome(), mConfig.RequiresLiftover);
        return new DriverData(driver, purity, comparisonChromosome, checkTranscript);
    }
}
