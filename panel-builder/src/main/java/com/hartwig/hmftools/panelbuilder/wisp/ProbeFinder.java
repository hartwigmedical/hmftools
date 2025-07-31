package com.hartwig.hmftools.panelbuilder.wisp;

import static java.lang.String.format;

import static com.hartwig.hmftools.common.genome.gc.GcCalcs.calcGcPercent;
import static com.hartwig.hmftools.common.genome.refgenome.RefGenomeSource.loadRefGenome;
import static com.hartwig.hmftools.common.utils.config.ConfigUtils.addLoggingOptions;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.wisp.common.CommonUtils.APP_NAME;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_CATEGORY;
import static com.hartwig.hmftools.wisp.common.CommonUtils.FLD_VARIANT;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;

import com.google.common.collect.Lists;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;
import com.hartwig.hmftools.common.utils.config.ConfigBuilder;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.jetbrains.annotations.NotNull;

public class ProbeFinder
{
    private final ProbeConfig mConfig;
    private final List<Variant> mCommonVariants;
    private final RefGenomeInterface mRefGenome;
    private final BufferedWriter mWriter;

    private static final Logger LOGGER = LogManager.getLogger(ProbeFinder.class);

    public ProbeFinder(final ConfigBuilder configBuilder)
    {
        mConfig = new ProbeConfig(configBuilder);
        mRefGenome = loadRefGenome(mConfig.RefGenomeFile);

        mCommonVariants = Lists.newArrayList();

        mWriter = initialiseWriter();
    }

    public void run()
    {
        if(mRefGenome == null)
        {
            System.exit(1);
        }

        LOGGER.info("sample({}) running probe variant selection", mConfig.SampleId);

        if(mConfig.ReferenceVariantsFile != null)
        {
            List<Variant> referenceVariants = ReferenceMutation.loadKnownMutations(mConfig.ReferenceVariantsFile);

            if(referenceVariants == null)
            {
                System.exit(1);
            }

            mCommonVariants.addAll(referenceVariants);
        }

        processSample(mConfig.SampleId);

        closeBufferedWriter(mWriter);

        LOGGER.info("Probe variation selection complete");
    }

    private void processSample(String sampleId)
    {
        List<Variant> selectedVariants = selectSampleVariants(sampleId);
        writeVariants(selectedVariants);
    }

    private List<Variant> selectSampleVariants(final String sampleId)
    {
        mConfig.checkSampleDirectories(sampleId);

        try
        {
            List<Variant> variants = Lists.newArrayList();
            variants.addAll(mCommonVariants);
            variants.addAll(SomaticMutation.loadSomatics(sampleId, mConfig));
            variants.addAll(GermlineMutation.loadGermlineMutations(sampleId, mConfig));
            variants.addAll(StructuralVariant.loadStructuralVariants(sampleId, mConfig));
            variants.addAll(GermlineSv.loadGermlineStructuralVariants(sampleId, mConfig));

            variants.forEach(x -> x.generateSequences(mRefGenome, mConfig));

            return VariantSelection.selectVariants(variants, mConfig);
        }
        catch(Exception e)
        {
            String error = format("sample data loading error: %s", e);
            LOGGER.error(e);
            throw new RuntimeException(error);
        }
    }

    private BufferedWriter initialiseWriter()
    {
        try
        {
            String filename = mConfig.OutputDir;

            filename += mConfig.SampleId + ".probe_variants";

            if(mConfig.OutputId != null)
            {
                filename += "." + mConfig.OutputId;
            }

            filename += TSV_EXTENSION;

            BufferedWriter writer = createBufferedWriter(filename, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            sj.add(CommonUtils.FLD_CATEGORY)
                    .add("Status")
                    .add(CommonUtils.FLD_VARIANT)
                    .add("Reported")
                    .add("CopyNumber")
                    .add("Vaf")
                    .add("TumorFrags");
            sj.add("PhasedVariants").add("Gene").add("Type").add("Sequence").add("GcPercent").add("OtherData");

            writer.write(sj.toString());
            writer.newLine();
            return writer;
        }
        catch(IOException e)
        {
            LOGGER.error("failed to initialise output file: {}", e.toString());
            return null;
        }
    }

    private void writeVariants(final List<Variant> selectedVariants)
    {
        try
        {
            for(Variant variant : selectedVariants)
            {
                if(!mConfig.WriteAll && !variant.isSelected())
                {
                    continue;
                }

                StringJoiner variantInfo = new StringJoiner(TSV_DELIM);

                variantInfo.add(variant.categoryType().toString());
                variantInfo.add(variant.selectionStatus().toString());
                variantInfo.add(variant.description());
                variantInfo.add(String.valueOf(variant.reported()));
                variantInfo.add(format("%.2f", variant.copyNumber()));
                variantInfo.add(format("%.2f", variant.vaf()));
                variantInfo.add(String.valueOf(variant.tumorFragments()));
                variantInfo.add(String.valueOf(variant.hasPhaseVariants()));
                variantInfo.add(variant.gene());

                StringJoiner sj = new StringJoiner(TSV_DELIM);
                sj.add(variantInfo.toString());
                sj.add("ALT");
                sj.add(variant.sequence());
                sj.add(format("%.2f", variant.gc()));
                sj.add(variant.otherData());

                mWriter.write(sj.toString());
                mWriter.newLine();

                for(String refSequence : variant.refSequences())
                {
                    StringJoiner refSj = new StringJoiner(TSV_DELIM);
                    refSj.add(variantInfo.toString());
                    refSj.add("REF");
                    refSj.add(refSequence);
                    refSj.add(format("%.2f", calcGcPercent(refSequence)));
                    refSj.add("");
                    mWriter.write(refSj.toString());
                    mWriter.newLine();
                }
            }
        }
        catch(IOException e)
        {
            LOGGER.error(" failed to write variants: {}", e.toString());
        }
    }

    public static void main(@NotNull final String[] args)
    {
        ConfigBuilder configBuilder = new ConfigBuilder(CommonUtils.APP_NAME);
        ProbeConfig.addConfig(configBuilder);

        addLoggingOptions(configBuilder);

        configBuilder.checkAndParseCommandLine(args);

        ProbeFinder application = new ProbeFinder(configBuilder);
        application.run();
    }
}
