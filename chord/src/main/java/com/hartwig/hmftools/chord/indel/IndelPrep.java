package com.hartwig.hmftools.chord.indel;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.prep.MutTypeCount;
import com.hartwig.hmftools.chord.variant.SmallVariant;
import com.hartwig.hmftools.chord.variant.VcfFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class IndelPrep
{
    private final ChordConfig mConfig;

    @NotNull private final RefGenomeSource mRefGenome;

    List<IndelDetails> mIndelDetailsList = new ArrayList<>();

    private static final String INDEL_DETAILS_FILE_SUFFIX = ".chord.indel.details.tsv";

    public IndelPrep(ChordConfig config)
    {
        mConfig = config;

        mRefGenome = RefGenomeSource.loadRefGenome(config.RefGenomeFile);
    }

    private List<SmallVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.purpleSomaticVcfFile(sampleId), mConfig.IncludeNonPass);
        List<VariantContext> variantContexts = vcfFile.loadVariants();

        List<SmallVariant> variants = new ArrayList<>();
        for(VariantContext variantContext : variantContexts)
        {
            SmallVariant smallVariant = new SmallVariant(variantContext);

            if(!smallVariant.isIndel())
                continue;

            variants.add(smallVariant);
        }

        return variants;
    }

    public List<MutTypeCount> extractSampleData(String sampleId)
    {
        try
        {
            CHORD_LOGGER.info("Counting indel contexts from sample: {}", sampleId);

            List<SmallVariant> variants = loadVariants(sampleId);

            Map<String, Integer> contextCountsMap = IndelContext.initializeCounts();
            for(SmallVariant variant : variants)
            {
                IndelVariant indel = new IndelVariant(variant);

                IndelDetails indelDetails = IndelDetails.from(indel, mRefGenome);
                if(mConfig.WriteDetailedFiles)
                {
                    mIndelDetailsList.add(indelDetails);
                }

                IndelContext indelContext = IndelContext.from(indelDetails);
                String indelContextName = indelContext.getContextName();

                contextCountsMap.compute(indelContextName, (k,v) -> v + 1);
            }

            if(mConfig.WriteDetailedFiles)
            {
                String indelDetailsPath = mConfig.OutputDir + "/" + sampleId + INDEL_DETAILS_FILE_SUFFIX;
                CHORD_LOGGER.info("Writing indel details to: {}", indelDetailsPath);
                writeDetails(indelDetailsPath, mIndelDetailsList);
            }

            List<MutTypeCount> contextCountsList = new ArrayList<>();
            for(String indelContextName : contextCountsMap.keySet())
            {
                int count = contextCountsMap.get(indelContextName);
                MutTypeCount mutTypeCount = new MutTypeCount(indelContextName, count);

                CHORD_LOGGER.debug(mutTypeCount);

                contextCountsList.add(mutTypeCount);
            }

            return contextCountsList;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("sample({}) failed to count indel contexts:", sampleId);
            e.printStackTrace();
            System.exit(1);
            return null;
        }
    }

    private static void writeDetails(String path, List<IndelDetails> indelDetailsList) throws IOException
    {
        BufferedWriter writer = IndelDetails.initializeWriter(path);

        for(IndelDetails indelDetails : indelDetailsList)
            indelDetails.writeLine(writer);

        writer.close();
    }

}
