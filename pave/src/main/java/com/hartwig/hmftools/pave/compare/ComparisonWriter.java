package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.utils.file.FileDelimiters.TSV_EXTENSION;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.impact.CodingContext;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.impact.ProteinContext;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.impact.VariantTransImpact;

public class ComparisonWriter
{
    private final BufferedWriter mImpactWriter;
    private final BufferedWriter mTransImpactWriter;
    private final BufferedWriter mRefVariantWriter;
    private final GeneDataCache mGeneDataCache;

    public ComparisonWriter(final GeneDataCache geneDataCache, final ComparisonConfig config)
    {
        mGeneDataCache = geneDataCache;
        mImpactWriter = initialiseImpactWriter(config.OutputDir, config.OutputId);
        mTransImpactWriter = config.WriteTransData ? initialiseTransWriter(config.OutputDir, config.OutputId) : null;

        mRefVariantWriter = config.ReferenceVariantsFile == null ? initialiseRefVariantWriter(config.OutputDir, config.OutputId) : null;
    }

    public void close()
    {
        closeBufferedWriter(mImpactWriter);
        closeBufferedWriter(mTransImpactWriter);
        closeBufferedWriter(mRefVariantWriter);
    }

    private static String formFilename(final String outputDir, final String outputId, final String fileType)
    {
        if(outputId != null)
            return outputDir + "pave_compare_" + fileType + "_" +  outputId + TSV_EXTENSION;
        else
            return outputDir + "pave_compare_" + fileType + TSV_EXTENSION;
    }

    public synchronized void writeVariantData(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant, boolean hasDiff)
    {
        writeImpactData(sampleId, variant, variantImpact, refVariant);
        writeTransImpactData(sampleId, variant, refVariant);

        if(hasDiff)
            writeRefVariant(sampleId, refVariant);
    }

    private BufferedWriter initialiseImpactWriter(final String outputDir, final String outputId)
    {
        try
        {
            String fileName = formFilename(outputDir, outputId, "impacts");

            PV_LOGGER.info("writing impact comparison file: {}", fileName);

            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId\t");
            writer.write(VariantData.tsvHeader());
            writer.write("\tGeneName\tIsDriver\tCanonEffects\tCanonCodingEffect\tHgvsCoding\tHgvsProtein");
            writer.write("\tWorstCodingEffect\tReported\tGenesAffected");
            writer.write("\tOrigGeneName\tOrigCanonEffects\tOrigCanonCodingEffect");
            writer.write("\tOrigWorstCodingEffect\tOrigHgvsCoding\tOrigHgvsProtein\tOrigReported");
            writer.newLine();

            return writer;
        }
        catch(IOException e)
        {
            PV_LOGGER.error("failed to initialise CSV file output: {}", e.toString());
            return null;
        }
    }

    private void writeImpactData(
            final String sampleId, final VariantData variant, final VariantImpact variantImpact, final RefVariantData refVariant)
    {
        if(mImpactWriter == null)
            return;

        try
        {
            mImpactWriter.write(String.format("%s\t%s\t%s",
                    sampleId, variant.tsvData(), variantImpact.CanonicalGeneName));

            boolean isDriver = mGeneDataCache.isDriverPanelGene(variantImpact.CanonicalGeneName);

            mImpactWriter.write(String.format("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d",
                    isDriver, variantImpact.CanonicalEffect, variantImpact.CanonicalCodingEffect,
                    variantImpact.CanonicalHgvsCoding, variantImpact.CanonicalHgvsProtein,
                    variantImpact.WorstCodingEffect, variant.reported(), variantImpact.GenesAffected));

            mImpactWriter.write(String.format("\t%s\t%s\t%s\t%s\t%s\t%s\t%s",
                    refVariant.Gene, refVariant.CanonicalEffect, refVariant.CanonicalCodingEffect, refVariant.WorstCodingEffect,
                    refVariant.HgvsCodingImpact, refVariant.HgvsProteinImpact, refVariant.Reported));

            mImpactWriter.newLine();
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

}
