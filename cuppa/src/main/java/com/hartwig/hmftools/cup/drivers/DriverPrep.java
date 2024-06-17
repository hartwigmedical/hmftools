package com.hartwig.hmftools.cup.drivers;

import static com.hartwig.hmftools.cup.prep.CategoryType.DRIVER;
import static com.hartwig.hmftools.common.variant.VariantType.INDEL;
import static com.hartwig.hmftools.common.variant.msi.MicrosatelliteStatus.MSS;
import static com.hartwig.hmftools.common.virus.VirusLikelihoodType.HIGH;
import static com.hartwig.hmftools.common.virus.VirusLikelihoodType.UNKNOWN;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.DNA;

import java.io.File;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;

import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalog;
import com.hartwig.hmftools.common.drivercatalog.DriverCatalogFile;
import com.hartwig.hmftools.common.drivercatalog.DriverType;
import com.hartwig.hmftools.common.fusion.KnownFusionType;
import com.hartwig.hmftools.common.linx.FusionLikelihoodType;
import com.hartwig.hmftools.common.linx.LinxFusion;
import com.hartwig.hmftools.common.purple.PurityContext;
import com.hartwig.hmftools.common.purple.PurityContextFile;
import com.hartwig.hmftools.common.virus.AnnotatedVirus;
import com.hartwig.hmftools.common.virus.AnnotatedVirusFile;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;
import com.hartwig.hmftools.cup.somatics.SomaticVariantsLoader;
import com.hartwig.hmftools.cup.somatics.SomaticVariant;

public class DriverPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    LinkedHashMap<DataItem.Index, DataItem> mDataItemsMap = new LinkedHashMap<>();

    private static final double DRIVER_PRESENT_LIKELIHOOD = 1.0;
    private static final String FLOAT_FORMAT_LIKELIHOOD = "%.4f";

    private static final String AMP_SUFFIX = ".amp";
    private static final String MUTATION_SUFFIX = ".mut";
    private static final String INDEL_SUFFIX = ".indel";

    private static final String PROMISCUOUS_5_SUFFIX = "_PROM5";
    private static final String PROMISCUOUS_3_SUFFIX = "_PROM3";

    private static final int INDEL_MAX_REPEAT_COUNT = 6;
    private static final String INDEL_ALB = "ALB";
    private static final String INDEL_SFTPB = "SFTPB";
    private static final String INDEL_SLC34A2 = "SLC34A2";

    public DriverPrep(final PrepConfig config)
    {
        mConfig = config;
    }

    @Override
    public CategoryType categoryType() { return DRIVER; }

    private void addDataItem(DataItem dataItem)
    {
        if(!mDataItemsMap.containsKey(dataItem.Index))
        {
            mDataItemsMap.put(dataItem.Index, dataItem);
        }

        // De-duplicate features by max likelihood
        float newDataItemValue = Float.parseFloat(dataItem.Value);
        float existingDataItemValue = Float.parseFloat(mDataItemsMap.get(dataItem.Index).Value);

        if(newDataItemValue > existingDataItemValue)
        {
            mDataItemsMap.put(dataItem.Index, dataItem);
        }
    }

    public String getDriverCatalogFile(final String sampleId) throws NoSuchFileException
    {
        String path;

        path = mConfig.linxDriverCatalogFile(sampleId);
        if(new File(path).isFile())
            return path;

        path = mConfig.purpleDriverCatalogFile(sampleId);
        if(new File(path).isFile())
            return path;

        CUP_LOGGER.error("sample({}) has no LINX or PURPLE driver catalog file", sampleId);
        throw new NoSuchFileException(path);
    }

    public void getDriverMutationsFromCatalog(String sampleId) throws IOException
    {
        String driverCatalogFile = getDriverCatalogFile(sampleId);
        final List<DriverCatalog> drivers = DriverCatalogFile.read(driverCatalogFile);

        for(final DriverCatalog driver : drivers)
        {
            if(DriverType.isGermline(driver.driver()))
                continue;

            double likelihood = driver.driverLikelihood();

            String featureName = driver.gene();
            if(driver.driver() == DriverType.AMP || driver.driver() == DriverType.PARTIAL_AMP)
            {
                featureName += AMP_SUFFIX;
            } else {
                featureName += MUTATION_SUFFIX;
            }

            DataItem dataItem = new DataItem(DNA, ItemType.DRIVER, featureName, likelihood, FLOAT_FORMAT_LIKELIHOOD);
            addDataItem(dataItem);
        }
    }

    public void getRepeatIndelDrivers(String sampleId) throws IOException
    {
        PurityContext purityContext = PurityContextFile.readWithQC(
                mConfig.purpleQcFile(sampleId),
                mConfig.purplePurityFile(sampleId)
        );

        boolean isMicrosatelliteStable = purityContext.microsatelliteStatus() == MSS;

        final List<SomaticVariant> variants = SomaticVariantsLoader.loadFromConfig(mConfig, sampleId, null);
        for(SomaticVariant variant : variants)
        {
            String gene = variant.Gene;

            // TODO: Look for known hotspot mutations in CupConstants. Potentially rename method
            // KNOWN_MUTATIONS.stream().anyMatch(x -> x.matches(gene, variant.Type, variant.Ref, variant.Alt, variant.Position));

            boolean isKnownIndelGene = gene.equals(INDEL_ALB) || gene.equals(INDEL_SFTPB) || gene.equals(INDEL_SLC34A2);
            boolean isRepeatIndelDriver = isKnownIndelGene && variant.Type == INDEL && variant.RepeatCount <= INDEL_MAX_REPEAT_COUNT;

            if(isMicrosatelliteStable && isRepeatIndelDriver)
            {
                String featureName = gene + INDEL_SUFFIX;
                DataItem dataItem = new DataItem(DNA, ItemType.DRIVER, featureName, DRIVER_PRESENT_LIKELIHOOD, FLOAT_FORMAT_LIKELIHOOD);
                addDataItem(dataItem);
            }
        }
    }

    public void getFusions(String sampleId) throws IOException
    {
        final String fusionsFilename = mConfig.linxFusionFile(sampleId);
        List<LinxFusion> fusions = LinxFusion.read(fusionsFilename);

        for(final LinxFusion fusion : fusions)
        {
            if(!fusion.reported())
                continue;

            boolean isPromiscuous5 = fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_5.toString());
            boolean isPromiscuous3 = fusion.reportedType().equals(KnownFusionType.PROMISCUOUS_3.toString());

            final String fusionName;

            //TODO: split PROMISCUOUS_BOTH fusions into 2 sepearate DataItem entries?
            if(isPromiscuous5 || isPromiscuous3){

                if(fusion.likelihood() != FusionLikelihoodType.HIGH)
                    continue;

                final String[] genes = fusion.name().split("_");

                fusionName = isPromiscuous5 ?
                        genes[0] + PROMISCUOUS_5_SUFFIX :
                        genes[1] + PROMISCUOUS_3_SUFFIX;
            }
            else
            {
                fusionName = fusion.name();
            }

            DataItem dataItem = new DataItem(DNA, ItemType.FUSION, fusionName, DRIVER_PRESENT_LIKELIHOOD, FLOAT_FORMAT_LIKELIHOOD);
            addDataItem(dataItem);
        }
    }

    public void getVirusAnnotations(String sampleId) throws IOException
    {
        String viralAnnotationFilename = mConfig.viralAnnotationFile(sampleId);
        final List<AnnotatedVirus> virusAnnotations = new ArrayList<>();

        AnnotatedVirusFile.read(viralAnnotationFilename).stream().filter(AnnotatedVirus::reported).forEach(virusAnnotations::add);

        if(virusAnnotations.size() == 0)
            return;

        for(AnnotatedVirus annotatedVirus : virusAnnotations)
        {
            if(!annotatedVirus.reported())
                continue;

            // `virusDriverLikelihoodType() == UNKNOWN` does not mean that the likelihood is unknown, but that the virus is not a known
            // carcinogenic virus (i.e. annotatedVirus.reported()==false)
            double likelihood = (annotatedVirus.virusDriverLikelihoodType() == HIGH || annotatedVirus.virusDriverLikelihoodType() == UNKNOWN) ? 1 : 0.5;
            String virusName = ViralInsertionType.fromVirusName(annotatedVirus.name()).toString();

            DataItem dataItem = new DataItem(DNA, ItemType.VIRUS, virusName, likelihood, FLOAT_FORMAT_LIKELIHOOD);
            addDataItem(dataItem);
        }
    }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        try {
            getDriverMutationsFromCatalog(sampleId);
            getRepeatIndelDrivers(sampleId);
            getFusions(sampleId);
            getVirusAnnotations(sampleId);
        }
        catch(Exception e)
        {
            CUP_LOGGER.error("sample({}) failed to extract category({}):", sampleId, categoryType());
            e.printStackTrace();
            System.exit(1);
        }

        return new ArrayList<>(mDataItemsMap.values());
    }
}
