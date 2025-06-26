package com.hartwig.hmftools.chord.indel;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.prep.LoggingOptions;
import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.prep.SmallVariant;
import com.hartwig.hmftools.chord.prep.VariantTypePrep;
import com.hartwig.hmftools.chord.prep.VcfFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;

public class IndelPrep implements VariantTypePrep<IndelVariant>, LoggingOptions
{
    private final ChordConfig mConfig;
    private String mLogPrefix = "";

    private final RefGenomeSource mRefGenome;

    private final List<IndelDetails> mIndelDetailsList = new ArrayList<>();
    private final Map<Type, Integer> mSkippedVariantTypeCounts = new HashMap<>();


    private static final String INDEL_DETAILS_FILE_SUFFIX = ".chord.indel.details.tsv";

    public IndelPrep(ChordConfig config)
    {
        mConfig = config;

        if(config.RefGenomeFile == null)
            throw new IllegalStateException("Ref genome must be provided to run " + this.getClass().getSimpleName());

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
        List<VariantContext> variants = vcfFile.loadVariants();

        List<IndelVariant> indels = new ArrayList<>();

        for(VariantContext variantContext : variants)
        {
            if(!variantContext.isIndel())
            {
                if(CHORD_LOGGER.isDebugEnabled())
                {
                    mSkippedVariantTypeCounts.put(
                            variantContext.getType(),
                            mSkippedVariantTypeCounts.getOrDefault(variantContext.getType(), 0) + 1
                    );
                }

                CHORD_LOGGER.trace("{}Skipped variant: {}", mLogPrefix, variantContext);

                continue;
            }

            SmallVariant smallVariant = new SmallVariant(variantContext);
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

            if(!mSkippedVariantTypeCounts.isEmpty())
            {
                CHORD_LOGGER.debug("{}Skipped variants: {}", mLogPrefix, mSkippedVariantTypeCounts);
            }

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
