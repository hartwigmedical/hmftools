package com.hartwig.hmftools.chord.snv;

import static com.hartwig.hmftools.chord.ChordConstants.CHORD_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.NoSuchFileException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.prep.LoggingOptions;
import com.hartwig.hmftools.chord.prep.MutContextCount;
import com.hartwig.hmftools.chord.prep.SmallVariant;
import com.hartwig.hmftools.chord.prep.VariantTypePrep;
import com.hartwig.hmftools.chord.prep.VcfFile;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource;
import com.hartwig.hmftools.common.sigs.SnvSigUtils;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContext.Type;

public class SnvPrep implements VariantTypePrep<SmallVariant>, LoggingOptions
{
    private final ChordConfig mConfig;
    private String mLogPrefix = "";

    private final RefGenomeSource mRefGenome;

    private final List<SnvDetails> mSnvDetailsList = new ArrayList<>();
    private final Map<Type, Integer> mSkippedVariantTypeCounts = new HashMap<>();

    private static final String SNV_DETAILS_FILE_SUFFIX = ".chord.snv.details.tsv";

    private static final Set<Character> VALID_NUCLEOTIDES = Set.of(
            'A', 'C', 'G', 'T',
            'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N');

    private static final Set<Character> AMBIGUOUS_NUCLEOTIDES = Set.of(
            'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', 'N');

    private static final float MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION = 0.1f;

    public SnvPrep(ChordConfig config)
    {
        mConfig = config;

        if(config.RefGenomeFile == null)
            throw new IllegalStateException("Ref genome must be provided to run " + this.getClass().getSimpleName());

        mRefGenome = RefGenomeSource.loadRefGenome(config.RefGenomeFile);
    }

    @Override
    public SnvPrep logPrefix(String logPrefix)
    {
        mLogPrefix = logPrefix;
        return this;
    }

    @Override
    public List<SmallVariant> loadVariants(String sampleId) throws NoSuchFileException
    {
        VcfFile vcfFile = new VcfFile(mConfig.snvIndelVcfFile(sampleId), mConfig.IncludeNonPass).logPrefix(mLogPrefix);

        List<VariantContext> variants = vcfFile.loadVariants();

        List<SmallVariant> snvs = new ArrayList<>();
        for(VariantContext variantContext : variants)
        {
            if(!variantContext.isSNP())
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
    public List<MutContextCount> countMutationContexts(String sampleId)
    {
        try
        {
            CHORD_LOGGER.debug("{}Running {} - counting trinucleotide contexts", mLogPrefix, this.getClass().getSimpleName());

            List<SmallVariant> snvs = loadVariants(sampleId);
            CHORD_LOGGER.debug("{}Found {} SNVs", mLogPrefix, snvs.size());

            if(!mSkippedVariantTypeCounts.isEmpty())
            {
                CHORD_LOGGER.debug("{}Skipped variants: {}", mLogPrefix, mSkippedVariantTypeCounts);
            }

            Map<String,Integer> triNucNameCountsMap = initializeCounts();
            int ambiguousTriNucContextCount = 0;

            for(SmallVariant snv : snvs)
            {
                SnvDetails snvDetails = SnvDetails.from(sampleId, snv, mRefGenome);

                if(mConfig.WriteDetailedFiles)
                {
                    mSnvDetailsList.add(snvDetails);
                }

                if(triNucNameCountsMap.containsKey(snvDetails.mTriNucContext))
                {
                    triNucNameCountsMap.compute(snvDetails.mTriNucContext, (k, v) -> v + 1);
                }
                else if(isAmbiguousTriNucContext(snvDetails))
                {
                    CHORD_LOGGER.warn("{}Found ambiguous nucleotide in trinucleotide context '{}' for variant '{}:{}:{}:{}'",
                            mLogPrefix, snvDetails.mTriNucContext,
                            snvDetails.mChromosome, snvDetails.mPosition, snvDetails.mRefBases, snv.AltBases
                    );
                    ambiguousTriNucContextCount += 1;
                }
                else
                {
                    CHORD_LOGGER.error("{}Found invalid trinucleotide context '{}' for variant '{}:{}:{}:{}'",
                            mLogPrefix, snvDetails.mTriNucContext,
                            snvDetails.mChromosome, snvDetails.mPosition, snvDetails.mRefBases, snv.AltBases
                    );
                    return null;
                }
            }

            if(mConfig.WriteDetailedFiles)
            {
                String snvDetailsPath = mConfig.OutputDir + "/" + sampleId + SNV_DETAILS_FILE_SUFFIX;
                CHORD_LOGGER.debug("{}Writing SNV details to: {}", mLogPrefix, snvDetailsPath);
                writeDetails(snvDetailsPath, mSnvDetailsList);
            }

            if(ambiguousTriNucContextCount > MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION * snvs.size())
            {
                CHORD_LOGGER.error("{}Fraction of trinucleotide contexts with ambiguous bases '{}' exceeds maximum '{}'",
                        mLogPrefix,
                        String.format("%.2f", ambiguousTriNucContextCount / (float) snvs.size()),
                        String.format("%.2f", MAX_AMBIGUOUS_TRI_NUCLEOTIDE_CONTEXT_FRACTION)
                );
                return null;
            }

            List<MutContextCount> triNucCountsList = new ArrayList<>();

            for(String binName : triNucNameCountsMap.keySet())
            {
                String newBinName = SnvDetails.renameTriNucBin(binName);
                int count = triNucNameCountsMap.get(binName);

                MutContextCount mutTypeCount = new MutContextCount(newBinName, count);

                CHORD_LOGGER.trace(mutTypeCount);

                triNucCountsList.add(mutTypeCount);
            }
            CHORD_LOGGER.debug("{}Completed {}", mLogPrefix, this.getClass().getSimpleName());

            return triNucCountsList;
        }
        catch(Exception e)
        {
            CHORD_LOGGER.error("{}{} failed: {}", mLogPrefix, this.getClass().getSimpleName(), e.toString());
            e.printStackTrace();
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

    private boolean isAmbiguousTriNucContext(SnvDetails snvDetails)
    {
        String triNucSequence = snvDetails.mTriNucSequence;

        List<Character> relevantBases = List.of(
                snvDetails.mAltBases.charAt(0),
                triNucSequence.charAt(0),
                triNucSequence.charAt(1),
                triNucSequence.charAt(2)
        );

        return VALID_NUCLEOTIDES.containsAll(relevantBases) && relevantBases.stream().anyMatch(AMBIGUOUS_NUCLEOTIDES::contains);
    }
}