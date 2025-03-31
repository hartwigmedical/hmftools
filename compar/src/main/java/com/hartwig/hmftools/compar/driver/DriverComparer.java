package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.driver.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.compar.common.Category.DRIVER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.common.CommonUtils.determineComparisonChromosome;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKE_METHOD;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_MAX_COPY_NUMBER;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_MIN_COPY_NUMBER;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.driver.DriverCatalog;
import com.hartwig.hmftools.common.driver.DriverCatalogFile;
import com.hartwig.hmftools.common.linx.LinxDriver;
import com.hartwig.hmftools.compar.common.Category;
import com.hartwig.hmftools.compar.common.CommonUtils;
import com.hartwig.hmftools.compar.ComparConfig;
import com.hartwig.hmftools.compar.ComparableItem;
import com.hartwig.hmftools.compar.common.DiffThresholds;
import com.hartwig.hmftools.compar.common.FileSources;
import com.hartwig.hmftools.compar.common.Mismatch;
import com.hartwig.hmftools.compar.ItemComparer;
import com.hartwig.hmftools.patientdb.dao.DatabaseAccess;

public class DriverComparer implements ItemComparer
{
    private final ComparConfig mConfig;

    public DriverComparer(final ComparConfig config)
    {
        mConfig = config;
    }

    @Override
    public Category category() { return DRIVER; }

    @Override
    public void registerThresholds(final DiffThresholds thresholds)
    {
        thresholds.addFieldThreshold(FLD_LIKELIHOOD, 0.1, 0);
        thresholds.addFieldThreshold(FLD_MIN_COPY_NUMBER, 0.3, 0.15);
        thresholds.addFieldThreshold(FLD_MAX_COPY_NUMBER, 0.3, 0.15);
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_LIKE_METHOD, FLD_LIKELIHOOD, FLD_MIN_COPY_NUMBER, FLD_MAX_COPY_NUMBER);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);

        return drivers.stream().map(this::createDriverData).collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final String germlineSampleId, final FileSources fileSources)
    {
        final List<ComparableItem> comparableItems = Lists.newArrayList();

        try
        {
            // use Linx if present, otherwise Purple drivers
            String linxDriverFile = LinxDriver.generateCatalogFilename(fileSources.Linx, sampleId, true);
            String purpleDriverFile = DriverCatalogFile.generateSomaticFilename(fileSources.Purple, sampleId);

            List<DriverCatalog> drivers = Lists.newArrayList();

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
            String linxGermlineDriverFile = LinxDriver.generateCatalogFilename(fileSources.Linx, sampleId, false);

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

            for(DriverCatalog driver : drivers)
            {
                comparableItems.add(createDriverData(driver));
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load driver data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }

    private DriverData createDriverData(final DriverCatalog driver)
    {
        boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
        String comparisonChromosome = determineComparisonChromosome(driver.chromosome(), mConfig.RequiresLiftover);
        return new DriverData(driver, comparisonChromosome, checkTranscript);
    }
}
