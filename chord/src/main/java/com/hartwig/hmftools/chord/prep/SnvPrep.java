package com.hartwig.hmftools.chord.prep;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.variant.SmallVariant;
import com.hartwig.hmftools.chord.variant.VcfFile;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;
import com.hartwig.hmftools.common.variant.SageVcfTags;

import htsjdk.variant.variantcontext.VariantContext;

public class SnvPrep
{
    private final ChordConfig mConfig;

    public SnvPrep(ChordConfig config)
    {
        mConfig = config;
    }

    private List<VariantContext> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.purpleSomaticVcfFile(sampleId), mConfig.IncludeNonPass);
        return vcfFile.loadVariants();
    }

    private static Map<String,Integer> countTriNucContexts(List<VariantContext> variantContexts)
    {
        // Get bin names
        // Returns: C>A_ACA -> 0, C>A_ACC -> 1, C>A_ACG -> 2, etc. We only need the name (i.e. key)
        Map<String,Integer> triNucNameIndexMap = new LinkedHashMap<>();
        SnvSigUtils.populateBucketMap(triNucNameIndexMap);

        // Initialize count vector
        Map<String,Integer> triNucNameCountsMap = new LinkedHashMap<>();
        triNucNameIndexMap.keySet().forEach(i -> triNucNameCountsMap.put(i, 0));

        // Count trinucleotide contexts
        int snvCount = 0;
        for(VariantContext variantContext : variantContexts)
        {
            SmallVariant smallVariant = new SmallVariant(variantContext);

            if(!smallVariant.isSnv())
                continue;

            String triNucSequence = smallVariant.Context.getAttributeAsString(SageVcfTags.TRINUCLEOTIDE_CONTEXT, null);
            String triNucContext = SnvSigUtils.variantContext(smallVariant.RefBases, smallVariant.AltBases, triNucSequence);

            // CHORD_LOGGER.trace("{}:{}:{}:{} {}",
            //       variantContext.getContig(), variantContext.getStart(), refSeq, altSeq, renameTriNucBin(triNucContext));

            triNucNameCountsMap.compute(triNucContext, (k,v) -> v + 1);

            snvCount++;
        }

        CHORD_LOGGER.debug("Counted trinucleotide contexts for {} SNVs", snvCount);

        return triNucNameCountsMap;
    }

    private static String renameTriNucBin(String bucketName)
    {
        // Convert e.g. "C>A_ACA" to "A[C>A]A"
        // The latter is the format required by the CHORD random forest model

        String[] bucketNameSplit = bucketName.split("_");

        String substitutionType = bucketNameSplit[0];
        char[] triNucSequence = bucketNameSplit[1].toCharArray();

        return String.format("%s[%s]%s", triNucSequence[0], substitutionType, triNucSequence[2]);
    }

    private static List<MutTypeCount> makeCountsList(Map<String,Integer> triNucCounts)
    {
        List<MutTypeCount> counts = new ArrayList<>();

        for(String binName : triNucCounts.keySet())
        {
            String newBinName = renameTriNucBin(binName);
            int count = triNucCounts.get(binName);

            MutTypeCount mutTypeCount = new MutTypeCount(newBinName, count);

            CHORD_LOGGER.trace(mutTypeCount);

            counts.add(new MutTypeCount(newBinName, count));
        }

        return counts;
    }

    public List<MutTypeCount> extractSampleData(String sampleId)
    {
        try
        {
            CHORD_LOGGER.info("Counting SNV 96 trinucleotide contexts");

            List<VariantContext> variantContexts = loadVariants(sampleId);
            Map<String,Integer> triNucCountsMap = countTriNucContexts(variantContexts);
            List<MutTypeCount> triNucCountsList = makeCountsList(triNucCountsMap);

            return triNucCountsList;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("sample({}) failed to count SNV trinucleotide contexts:", sampleId);
            e.printStackTrace();
            System.exit(1);
            return null;
        }
    }
}
