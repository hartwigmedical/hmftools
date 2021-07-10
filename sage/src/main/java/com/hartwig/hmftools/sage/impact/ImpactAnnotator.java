package com.hartwig.hmftools.sage.impact;

import static com.hartwig.hmftools.common.ensemblcache.EnsemblDataCache.ENSEMBL_DATA_DIR;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.ConfigUtils.setLogLevel;
import static com.hartwig.hmftools.sage.SageCommon.SG_LOGGER;
import static com.hartwig.hmftools.sage.impact.ImpactConfig.REF_GENOME;

import static htsjdk.tribble.AbstractFeatureReader.getFeatureReader;

import java.io.IOException;
import java.util.List;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.ensemblcache.EnsemblGeneData;
import com.hartwig.hmftools.common.ensemblcache.TranscriptData;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.variant.SomaticVariant;
import com.hartwig.hmftools.common.variant.SomaticVariantFactory;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.jetbrains.annotations.NotNull;

import htsjdk.tribble.AbstractFeatureReader;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.filter.CompoundFilter;
import htsjdk.variant.vcf.VCFCodec;

public class ImpactAnnotator
{
    private final ImpactConfig mConfig;
    private final ImpactClassifier mImpactClassifier;
    private final GeneDataCache mGeneDataCache;

    public ImpactAnnotator(final CommandLine cmd)
    {
        mConfig = new ImpactConfig(cmd);

        mGeneDataCache = new GeneDataCache(cmd.getOptionValue(ENSEMBL_DATA_DIR), mConfig.RefGenVersion);

        RefGenomeInterface refGenome = loadRefGenome(cmd.getOptionValue(REF_GENOME));

        mImpactClassifier = new ImpactClassifier(refGenome);
    }

    public void run()
    {
        if(!mConfig.isValid())
        {
            SG_LOGGER.error("invalid config, exiting");
            System.exit(1);
        }

        SG_LOGGER.error("annotating variants with gene impacts");


        SG_LOGGER.error("annotation complete");
    }

    private void processVcfFile()
    {
        CompoundFilter filter = new CompoundFilter(true);
        SomaticVariantFactory variantFactory = new SomaticVariantFactory(filter);

        try
        {
            final AbstractFeatureReader<VariantContext, LineIterator> reader = getFeatureReader(
                    mConfig.VcfFile, new VCFCodec(), false);

            for (VariantContext variantContext : reader.iterator())
            {
                if (!filter.test(variantContext))
                    continue;

                // extract SnpEff data for comparison sake
                // SnpEffSummary
                // SnpEffAnnotation - per transcript

                final SomaticVariant variant = variantFactory.createVariant(mConfig.SampleId, variantContext).orElse(null);

                if(variant == null)
                    continue;

                processVariant(variant);
            }
        }
        catch(IOException e)
        {
            SG_LOGGER.error(" failed to read somatic VCF file({}): {}", mConfig.VcfFile, e.toString());
        }
    }

    private void processVariant(final SomaticVariant somaticVariant)
    {
        int variantPosition = (int)somaticVariant.position();
        List<EnsemblGeneData> geneCandidates = mGeneDataCache.findGenes(somaticVariant.chromosome(), variantPosition);

        if(geneCandidates.isEmpty())
        {
            // TODO: intergenic variant
            return;
        }

        // analyse against each of the genes and their transcripts
        for(EnsemblGeneData geneData : geneCandidates)
        {
            List<TranscriptData> transDataList = mGeneDataCache.findTranscripts(geneData.GeneId, variantPosition);

            transDataList.forEach(x -> mImpactClassifier.classifyVariant(x));
        }
    }

    public static void main(@NotNull final String[] args) throws ParseException
    {
        final Options options = ImpactConfig.createOptions();
        final CommandLine cmd = createCommandLine(args, options);

        setLogLevel(cmd);

        ImpactAnnotator impactAnnotator = new ImpactAnnotator(cmd);
        impactAnnotator.run();
    }

    @NotNull
    private static CommandLine createCommandLine(@NotNull final String[] args, @NotNull final Options options) throws ParseException
    {
        final CommandLineParser parser = new DefaultParser();
        return parser.parse(options, args);
    }

}
