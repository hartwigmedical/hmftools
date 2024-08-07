package com.hartwig.hmftools.compar.driver;

import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_GERMLINE;
import static com.hartwig.hmftools.common.drivercatalog.DriverType.DRIVERS_LINX_SOMATIC;
import static com.hartwig.hmftools.compar.common.Category.DRIVER;
import static com.hartwig.hmftools.compar.ComparConfig.CMP_LOGGER;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKELIHOOD;
import static com.hartwig.hmftools.compar.driver.DriverData.FLD_LIKE_METHOD;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

import com.google.common.collect.Lists;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverCategory;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.drivercatalog.ImmutableDriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.LikelihoodMethod;
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

import org.jetbrains.annotations.NotNull;

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
    }

    @Override
    public boolean processSample(final String sampleId, final List<Mismatch> mismatches)
    {
        return CommonUtils.processSample(this, mConfig, sampleId, mismatches);
    }

    @Override
    public List<String> comparedFieldNames()
    {
        return Lists.newArrayList(FLD_LIKE_METHOD, FLD_LIKELIHOOD);
    }

    @Override
    public List<ComparableItem> loadFromDb(final String sampleId, final DatabaseAccess dbAccess, final String sourceName)
    {
        final List<DriverCatalog> drivers = dbAccess.readDriverCatalog(sampleId);

        final List<ComparableItem> driverDataList = Lists.newArrayList();

        for(DriverCatalog driver : drivers)
        {
            boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
            driverDataList.add(new DriverData(driver, checkTranscript));
        }

        return driverDataList;
    }

    @Override
    public List<ComparableItem> loadFromOrangeJson(final String sampleId, final JsonObject json)
    {
        final List<ComparableItem> driverDataList = Lists.newArrayList();
        for(DriverCatalog driver : loadDriversFromOrangeJson(json))
        {
            boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
            driverDataList.add(new DriverData(driver, checkTranscript));
        }
        return driverDataList;
    }

    @NotNull
    private static List<DriverCatalog> loadDriversFromOrangeJson(final JsonObject json)
    {
        final List<DriverCatalog> drivers = Lists.newArrayList();
        for(JsonObject driverJson : extractPurpleDriversFromOrangeJson(json))
        {
            DriverCatalog driver = ImmutableDriverCatalog.builder()
                    .driverLikelihood(driverJson.get("driverLikelihood").getAsDouble())
                    .gene(driverJson.get("gene").getAsString())
                    .isCanonical(driverJson.get("isCanonical").getAsBoolean())
                    .likelihoodMethod(LikelihoodMethod.valueOf(driverJson.get("likelihoodMethod").getAsString()))
                    .transcript(driverJson.get("transcript").getAsString())
                    .driver(DriverType.valueOf(driverJson.get("type").getAsString()))
                    .chromosome("")  // unavailable in JSON
                    .chromosomeBand("")  // unavailable in JSON
                    .category(DriverCategory.ONCO)  // unavailable in JSON
                    .missense(-1)  // unavailable in JSON
                    .nonsense(-1)  // unavailable in JSON
                    .splice(-1)  // unavailable in JSON
                    .inframe(-1)  // unavailable in JSON
                    .frameshift(-1)  // unavailable in JSON
                    .biallelic(false)  // unavailable in JSON
                    .minCopyNumber(-1.)  // unavailable in JSON
                    .maxCopyNumber(-1.)  // unavailable in JSON
                    .build();
            drivers.add(driver);
        }

        for(JsonElement driverJson : json.getAsJsonObject("linx").getAsJsonArray("somaticDrivers"))
        {
            DriverType driverType = DriverType.checkConvertType(driverJson.getAsJsonObject().get("type").getAsString());
            if(DRIVERS_LINX_SOMATIC.contains(driverType))
            {
                DriverCatalog driver = ImmutableDriverCatalog.builder()
                        .gene(driverJson.getAsJsonObject().get("gene").getAsString())
                        .driver(driverType)
                        .chromosome("")  // unavailable in JSON
                        .chromosomeBand("")  // unavailable in JSON
                        .transcript("")  // unavailable in JSON
                        .isCanonical(true)  // unavailable in JSON
                        .category(DriverCategory.ONCO)  // unavailable in JSON
                        .likelihoodMethod(LikelihoodMethod.NONE)  // unavailable in JSON
                        .driverLikelihood(-1.)  // unavailable in JSON
                        .missense(-1)  // unavailable in JSON
                        .nonsense(-1)  // unavailable in JSON
                        .splice(-1)  // unavailable in JSON
                        .inframe(-1)  // unavailable in JSON
                        .frameshift(-1)  // unavailable in JSON
                        .biallelic(false)  // unavailable in JSON
                        .minCopyNumber(-1.)  // unavailable in JSON
                        .maxCopyNumber(-1.)  // unavailable in JSON
                        .build();
                drivers.add(driver);
            }
        }
        return drivers;
    }

    @NotNull
    private static List<JsonObject> extractPurpleDriversFromOrangeJson(final JsonObject json)
    {
        final JsonArray somaticDriverJsons = json.getAsJsonObject("purple").getAsJsonArray("somaticDrivers");
        final JsonArray germlineDriverJsons = json.getAsJsonObject("purple").getAsJsonArray("germlineDrivers");

        final JsonArray combinedDriverJsons = somaticDriverJsons.deepCopy();
        combinedDriverJsons.addAll(germlineDriverJsons);

        return StreamSupport.stream(combinedDriverJsons.spliterator(), true)
                .map(JsonElement::getAsJsonObject)
                .collect(Collectors.toList());
    }

    @Override
    public List<ComparableItem> loadFromFile(final String sampleId, final FileSources fileSources)
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
                boolean checkTranscript = mConfig.AlternateTranscriptDriverGenes.contains(driver.gene());
                comparableItems.add(new DriverData(driver, checkTranscript));
            }
        }
        catch(IOException e)
        {
            CMP_LOGGER.warn("sample({}) failed to load driver data: {}", sampleId, e.toString());
            return null;
        }

        return comparableItems;
    }
}
