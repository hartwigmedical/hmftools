package com.hartwig.hmftools.pave.compare;

import static com.hartwig.hmftools.common.utils.FileWriterUtils.closeBufferedWriter;
import static com.hartwig.hmftools.common.utils.FileWriterUtils.createBufferedWriter;
import static com.hartwig.hmftools.pave.PaveConfig.PV_LOGGER;

import java.io.BufferedWriter;
import java.io.IOException;

import com.hartwig.hmftools.common.variant.impact.VariantImpact;
import com.hartwig.hmftools.pave.CodingContext;
import com.hartwig.hmftools.pave.GeneDataCache;
import com.hartwig.hmftools.pave.ProteinContext;
import com.hartwig.hmftools.pave.VariantData;
import com.hartwig.hmftools.pave.VariantTransImpact;

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
            return outputDir + "pave_compare_" + fileType + "_" +  outputId + ".csv";
        else
            return outputDir + "pave_compare_" + fileType + ".csv";
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
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write("SampleId,");
            writer.write(VariantData.csvHeader());
            writer.write(",GeneName,IsDriver,CanonEffects,CanonCodingEffect,HgvsCoding,HgvsProtein");
            writer.write(",WorstCodingEffect,Reported,GenesAffected");
            writer.write(",SnpEffGeneName,SnpEffCanonEffects,SnpEffCanonCodingEffect");
            writer.write(",SnpEffWorstCodingEffect,SnpEffHgvsCoding,SnpEffHgvsProtein,SnpEffReported,SnpEffGenesAffected");
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
            mImpactWriter.write(String.format("%s,%s,%s",
                    sampleId, variant.csvData(), variantImpact.CanonicalGeneName));

            boolean isDriver = mGeneDataCache.isDriverPanelGene(variantImpact.CanonicalGeneName);

            mImpactWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%s,%d",
                    isDriver, variantImpact.CanonicalEffect, variantImpact.CanonicalCodingEffect,
                    variantImpact.CanonicalHgvsCoding, variantImpact.CanonicalHgvsProtein,
                    variantImpact.WorstCodingEffect, variant.reported(), variantImpact.GenesAffected));

            mImpactWriter.write(String.format(",%s,%s,%s,%s,%s,%s,%s,%d",
                    refVariant.Gene, refVariant.CanonicalEffect, refVariant.CanonicalCodingEffect, refVariant.WorstCodingEffect,
                    refVariant.HgvsCodingImpact, refVariant.HgvsProteinImpact, refVariant.Reported, refVariant.GenesAffected));

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
            BufferedWriter writer = createBufferedWriter(fileName, false);

            writer.write(String.format("SampleId,%s,GeneId,GeneName,%s,%s,%s",
                    VariantData.csvHeader(), VariantTransImpact.csvHeader(), CodingContext.csvHeader(), ProteinContext.csvHeader()));

            writer.write(",RealignHgvsCoding,RealignHgvsProtein");
            writer.write(",SnpEffCanonEffects,SnpEffHgvsCoding,SnpEffHgvsProtein,SagePhasedInframeIndel");

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
            // check for a name-mapping
            geneName = mGeneDataCache.getGeneNameFromSnpEff(refVariant.Gene);

            if(!variant.getImpacts().containsKey(geneName))
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

            mTransImpactWriter.write(String.format("%s,%s,%s,%s,%s,%s,%s",
                    sampleId, variant.csvData(), transImpact.TransData.GeneId, geneName,
                    transImpact.toCsv(), transImpact.codingContext().toCsv(),
                    transImpact.proteinContext() != null ? transImpact.proteinContext().toCsv() : ProteinContext.empty()));

            if(raImpact != null)
            {
                mTransImpactWriter.write(String.format(",%s,%s",
                        raImpact.codingContext().Hgvs, raImpact.proteinContext() != null ? raImpact.proteinContext().Hgvs : "N/A"));
            }
            else
            {
                mTransImpactWriter.write(",N/A,N/A");
            }

            mTransImpactWriter.write(String.format(",%s,%s,%s,%s",
                    refVariant.CanonicalEffect, refVariant.HgvsCodingImpact, refVariant.HgvsProteinImpact, refVariant.PhasedInframeIndel));

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
