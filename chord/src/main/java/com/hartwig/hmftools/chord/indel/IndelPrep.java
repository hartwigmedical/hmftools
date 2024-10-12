package com.hartwig.hmftools.chord.indel;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.common.LoggingOptions;
import com.hartwig.hmftools.chord.common.MutContextCount;
import com.hartwig.hmftools.chord.common.SmallVariant;
import com.hartwig.hmftools.chord.common.VariantTypePrep;
import com.hartwig.hmftools.chord.common.VcfFile;
import com.hartwig.hmftools.chord.snv.SnvPrep;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import org.jetbrains.annotations.NotNull;

import htsjdk.variant.variantcontext.VariantContext;

public class IndelPrep implements VariantTypePrep<IndelVariant>, LoggingOptions
{
    private final ChordConfig mConfig;
    private String mLogPrefix = "";

    @NotNull private final RefGenomeSource mRefGenome;

    List<IndelDetails> mIndelDetailsList = new ArrayList<>();

    private static final String INDEL_DETAILS_FILE_SUFFIX = ".chord.indel.details.tsv";

    public IndelPrep(ChordConfig config)
    {
        mConfig = config;

        if(config.RefGenomeFile == null)
            throw new IllegalArgumentException(this.getClass().getSimpleName() + " requires ref genome to be provided");

        mRefGenome = RefGenomeSource.loadRefGenome(config.RefGenomeFile);
    }

    @Override
    public IndelPrep logPrefix(String logPrefix)
    {
        mLogPrefix = logPrefix;
        return this;
    }

    @Override
    public List<IndelVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.snvIndelVcfFile(sampleId), mConfig.IncludeNonPass).logPrefix(mLogPrefix);
        List<VariantContext> variantContexts = vcfFile.loadVariants();

        List<IndelVariant> indels = new ArrayList<>();
        for(VariantContext variantContext : variantContexts)
        {
            SmallVariant smallVariant = new SmallVariant(variantContext);

            if(!smallVariant.isIndel())
                continue;

            IndelVariant indel = new IndelVariant(smallVariant);

            indels.add(indel);
        }

        return indels;
    }

    @Override
    public List<MutContextCount> countMutationContexts(String sampleId)
    {
        try
        {
            CHORD_LOGGER.debug("{}Running {} - counting microhomology and repeat contexts", mLogPrefix, this.getClass().getSimpleName());

            List<IndelVariant> indels = loadVariants(sampleId);
            CHORD_LOGGER.debug("{}Found {} indels", mLogPrefix, indels.size());

            Map<String, Integer> contextCountsMap = IndelContext.initializeCounts();

            for(IndelVariant indel : indels)
            {
                IndelDetails indelDetails = IndelDetails.from(sampleId, indel, mRefGenome);
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
                CHORD_LOGGER.debug("{}Writing indel details to: {}", mLogPrefix, indelDetailsPath);
                writeDetails(indelDetailsPath, mIndelDetailsList);
            }

            List<MutContextCount> contextCountsList = new ArrayList<>();
            for(String indelContextName : contextCountsMap.keySet())
            {
                int count = contextCountsMap.get(indelContextName);
                MutContextCount mutTypeCount = new MutContextCount(indelContextName, count);

                CHORD_LOGGER.trace(mutTypeCount);

                contextCountsList.add(mutTypeCount);
            }

            CHORD_LOGGER.debug("{}Completed {}", mLogPrefix, this.getClass().getSimpleName());

            return contextCountsList;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("{}{} failed: {}", mLogPrefix, this.getClass().getSimpleName(), e.toString());
            e.printStackTrace();
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
