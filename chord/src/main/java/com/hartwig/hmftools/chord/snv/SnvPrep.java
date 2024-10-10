package com.hartwig.hmftools.chord.snv;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.common.MutTypeCount;
import com.hartwig.hmftools.chord.common.SmallVariant;
import com.hartwig.hmftools.chord.common.VariantTypePrep;
import com.hartwig.hmftools.chord.common.VcfFile;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;

import htsjdk.variant.variantcontext.VariantContext;

public class SnvPrep implements VariantTypePrep<SmallVariant>
{
    private final ChordConfig mConfig;

    List<SnvDetails> mSnvDetailsList = new ArrayList<>();

    private static final String SNV_DETAILS_FILE_SUFFIX = ".chord.snv.details.tsv";

    public SnvPrep(ChordConfig config)
    {
        mConfig = config;
    }

    @Override
    public List<SmallVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.snvIndelVcfFile(sampleId), mConfig.IncludeNonPass);

        List<VariantContext> variants = vcfFile.loadVariants();

        List<SmallVariant> snvs = new ArrayList<>();
        for(VariantContext variantContext : variants)
        {
            SmallVariant smallVariant = new SmallVariant(variantContext);

            if(!smallVariant.isSnv())
                continue;

            snvs.add(smallVariant);
        }

        return snvs;
    }

    private static Map<String, Integer> initializeCounts()
    {
        // Get bin names
        // Returns: C>A_ACA -> 0, C>A_ACC -> 1, C>A_ACG -> 2, etc. We only need the name (i.e. key)
        Map<String,Integer> triNucNameIndexMap = new LinkedHashMap<>();
        SnvSigUtils.populateBucketMap(triNucNameIndexMap);

        // Initialize counts
        Map<String,Integer> triNucNameCountsMap = new LinkedHashMap<>();
        triNucNameIndexMap.keySet().forEach(i -> triNucNameCountsMap.put(i, 0));

        return triNucNameCountsMap;
    }

    @Override
    public List<MutTypeCount> countMutationContexts(String sampleId)
    {
        try
        {
            CHORD_LOGGER.info("Extracting SNV 96 trinucleotide contexts");

            List<SmallVariant> snvs = loadVariants(sampleId);
            CHORD_LOGGER.debug("Loaded {} SNVs", snvs.size());

            CHORD_LOGGER.debug("Initializing counts");
            Map<String,Integer> triNucNameCountsMap = initializeCounts();

            CHORD_LOGGER.debug("Populating counts");
            for(SmallVariant snv : snvs)
            {
                SnvDetails snvDetails = SnvDetails.from(sampleId, snv);

                if(mConfig.WriteDetailedFiles)
                    mSnvDetailsList.add(snvDetails);

                triNucNameCountsMap.compute(snvDetails.mTriNucContext, (k,v) -> v + 1);
            }

            if(mConfig.WriteDetailedFiles)
            {
                String snvDetailsPath = mConfig.OutputDir + "/" + sampleId + SNV_DETAILS_FILE_SUFFIX;
                CHORD_LOGGER.info("Writing SNV details to: {}", snvDetailsPath);
                writeDetails(snvDetailsPath, mSnvDetailsList);
            }

            List<MutTypeCount> triNucCountsList = new ArrayList<>();

            for(String binName : triNucNameCountsMap.keySet())
            {
                String newBinName = SnvDetails.renameTriNucBin(binName);
                int count = triNucNameCountsMap.get(binName);

                MutTypeCount mutTypeCount = new MutTypeCount(newBinName, count);

                CHORD_LOGGER.trace(mutTypeCount);

                triNucCountsList.add(mutTypeCount);
            }

            return triNucCountsList;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("sample({}) failed to count SNV trinucleotide contexts", sampleId);
            e.printStackTrace();
            System.exit(1);
            return null;
        }
    }

    private static void writeDetails(String path, List<SnvDetails> snvDetailsList) throws IOException
    {
        BufferedWriter writer = SnvDetails.initializeWriter(path);

        for(SnvDetails snvDetails : snvDetailsList)
            snvDetails.writeLine(writer);

        writer.close();
    }
}
