package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_DELIM;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.StringJoiner;
import java.util.stream.Collectors;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.VariantData;

public class ComparisonWriter
{
    private final ComparisonConfig mConfig;
    private final BufferedWriter mImpactDiffWriter;
    // private final BufferedWriter mTransImpactWriter;
    // private final BufferedWriter mRefVariantWriter;
    private final GeneDataCache mGeneDataCache;

    public ComparisonWriter(final GeneDataCache geneDataCache, final ComparisonConfig config)
    {
        mConfig = config;
        mGeneDataCache = geneDataCache;
        mImpactDiffWriter = initialiseDiffWriter();

        // mTransImpactWriter = config.WriteTransData ? initialiseTransWriter(config.OutputDir, config.OutputId) : null;
        /// mRefVariantWriter = config.ReferenceVariantsFile == null ? initialiseRefVariantWriter(config.OutputDir, config.OutputId) : null;
    }

    public void close()
    {
        closeBufferedWriter(mImpactDiffWriter);
        // closeBufferedWriter(mTransImpactWriter);
        // closeBufferedWriter(mRefVariantWriter);
    }

    private String formFilename(final String fileType)
    {
        String filename = mConfig.OutputDir + "pave_compare." + fileType;

        if(mConfig.OutputId != null)
            filename += "." +  mConfig.OutputId;

        filename += TSV_EXTENSION;
        return filename;
    }

    private BufferedWriter initialiseDiffWriter()
    {
        try
        {
            String fileName = formFilename("diffs");

            PV_LOGGER.info("writing impact comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(!mConfig.singleSample())
                sj.add("SampleId");

            sj.add(VariantData.tsvHeader());

            if(mConfig.LiftoverCache != null)
                sj.add("PositionV37");

            sj.add("GeneName").add("IsDriver").add("Diffs");
            sj.add("IsReportedRef").add("IsReportedNew");
            sj.add("CanonicalEffectsRef").add("CanonicalEffectsNew");
            sj.add("WorstCodingEffectRef").add("WorstCodingEffectNew");
            sj.add("CanonicalCodingEffectRef").add("CanonicalCodingEffectNew");
            sj.add("HgvsCodingRef").add("HgvsCodingNew");
            sj.add("HgvsProteinRef").add("HgvsProteinNew");

            writer.write(sj.toString());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise diff file output: {}", e.toString());
            return null;
        }
    }

    public synchronized void writeVariantDiff(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant,
            final boolean isDriverGene, final List<String> diffs)
    {
        if(mImpactDiffWriter == null)
            return;

        try
        {
            StringJoiner sj = new StringJoiner(TSV_DELIM);

            if(!mConfig.singleSample())
                sj.add(sampleId);

            sj.add(variant.tsvData());

            if(mConfig.LiftoverCache != null)
                sj.add(String.valueOf(refVariant.liftedPosition() != null ? refVariant.liftedPosition() : -1));

            String diffsStr = diffs.stream().collect(Collectors.joining(";"));

            sj.add(variantImpact.GeneName).add(String.valueOf(isDriverGene)).add(diffsStr);

            sj.add(String.valueOf(refVariant.Reported)).add(String.valueOf(variant.reported()));
            sj.add(refVariant.CanonicalEffect).add(variantImpact.CanonicalEffect);
            sj.add(String.valueOf(refVariant.WorstCodingEffect)).add(String.valueOf(variantImpact.WorstCodingEffect));
            sj.add(String.valueOf(refVariant.CanonicalCodingEffect)).add(String.valueOf(variantImpact.CanonicalCodingEffect));
            sj.add(refVariant.HgvsCodingImpact).add(variantImpact.CanonicalHgvsCoding);
            sj.add(refVariant.HgvsProteinImpact).add(variantImpact.CanonicalHgvsProtein);

            mImpactDiffWriter.write(sj.toString());
            mImpactDiffWriter.newLine();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise diff file output: {}", e.toString());
        }
    }

    /*
    private void writeImpactData(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        if(mImpactDiffWriter == null)
            return;

        try
        {
            mImpactDiffWriter.write(String.format("%s\t%s\t%s",
                    sampleId, variant.tsvData(), variantImpact.CanonicalGeneName));

            boolean isDriver = mGeneDataCache.isDriverPanelGene(variantImpact.CanonicalGeneName);

            mImpactDiffWriter.write(String.format("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d",
                    isDriver, variantImpact.CanonicalEffect, variantImpact.CanonicalCodingEffect,
                    variantImpact.CanonicalHgvsCoding, variantImpact.CanonicalHgvsProtein,
                    variantImpact.WorstCodingEffect, variant.reported(), variantImpact.GenesAffected));

            mImpactDiffWriter.write(String.format("\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    refVariant.Gene, refVariant.CanonicalEffect, refVariant.CanonicalCodingEffect, refVariant.WorstCodingEffect,
                    refVariant.HgvsCodingImpact, refVariant.HgvsProteinImpact, refVariant.Reported));

            mImpactDiffWriter.newLine();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write impact data: {}", e.toString());
            return;
        }
    }

    private BufferedWriter initialiseTransWriter(final String outputDir, final String outputId)
    {
        try
        {
            String fileName = formFilename(outputDir, outputId, "trans_impacts");

            PV_LOGGER.info("writing transcript comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(String.format("SampleId\t%s\tGeneId\tGeneName\t%s\t%s\t%s",
                    VariantData.tsvHeader(), VariantTransImpact.tsvHeader(), CodingContext.tsvHeader(), ProteinContext.tsvHeader()));

            writer.write("\tRealignHgvsCoding\tRealignHgvsProtein");
            writer.write("\tOrigCanonEffects\tOrigHgvsCoding\tOrigHgvsProtein");

            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise transcript impact file: {}", e.toString());
            return null;
        }
    }

    private void writeTransImpactData(final String sampleId, final VariantData variant, final RefVariantData refVariant)
    {
        if(mTransImpactWriter == null)
            return;

        String geneName = refVariant.Gene;

        if(!variant.getImpacts().containsKey(geneName))
        {
            if(!variant.getImpacts().containsKey(refVariant.Gene))
                return;
        }

        // only write if canonical transcript impact data can be found
        VariantTransImpact transImpact = variant.getImpacts().get(geneName).stream()
                .filter(x -> x.TransData.IsCanonical).findFirst().orElse(null);

        if(transImpact == null)
            return;

        try
        {
            VariantTransImpact raImpact = variant.getRealignedImpact(geneName, transImpact);

            mTransImpactWriter.write(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    sampleId, variant.tsvData(), transImpact.TransData.GeneId, geneName,
                    transImpact.toTsv(), transImpact.codingContext().toTsv(),
                    transImpact.proteinContext() != null ? transImpact.proteinContext().toTsv() : ProteinContext.empty()));

            if(raImpact != null)
            {
                mTransImpactWriter.write(String.format("\t%s\t%s",
                        raImpact.codingContext().Hgvs, raImpact.proteinContext() != null ? raImpact.proteinContext().Hgvs : "N/A"));
            }
            else
            {
                mTransImpactWriter.write("\tN/A\tN/A");
            }

            mTransImpactWriter.write(String.format("\t%s\t%s\t%s",
                    refVariant.CanonicalEffect, refVariant.HgvsCodingImpact, refVariant.HgvsProteinImpact));

            mTransImpactWriter.newLine();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write transcript impact data: {}", e.toString());
            return;
        }
    }

    private BufferedWriter initialiseRefVariantWriter(final String outputDir, final String outputId)
    {
        try
        {
            String fileName = outputId != null ?
                outputDir + "prod_ref_variants_" + outputId + ".tsv" : outputDir + "prod_ref_variants.tsv";

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(RefVariantData.tsvHeader());
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise ref variant file: {}", e.toString());
            return null;
        }
    }

    private void writeRefVariant(final String sampleId, final RefVariantData refVariant)
    {
        if(mRefVariantWriter == null)
            return;

        try
        {
            mRefVariantWriter.write(refVariant.tsvData(sampleId));
            mRefVariantWriter.newLine();
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to write ref variant data: {}", e.toString());
            return;
        }
    }
    */
}
